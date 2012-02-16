classdef RBJavaExport < handle
    % RBJavaExport: Exporting rbmatlab models for the JRB project.
    % This class is only here for having a safe version inside a git repo.
    %
    % Exports rbmatlab models to a specific folder.
    % 
    % The models can be exported in two formats, either for the JRB project (ans thus
    % ROMSim/Android) or for the rbappmit-Android app.
    %
    %
    % @author Daniel Wirtz @date 2011-08-26
    % Initial code for rbappmit-format export by Markus Dihlmann, 21.08.2011
    
    properties
        % The model export version to use. Can be either 0 for for JRB models or 1 for
        % rbappmit-models.
        %
        % @default 0
        Version = 0;
        
        % The folder to write the model data to.
        TargetFolder = '/home/dwirtz/aghwww/romsim/rbm_advec';
        
        % The source of the JRB software
        JRBSource = '/home/dwirtz/aghhome/Software/Eclipse/JRB/src';
        
        % Some additional export parameters.
        %
        % Optional fields of Params:
        %   model_type: defining 'system type' in input.in file. F.ex.
        %               'LINEAR_UNSTEADY'
        %   scm_type:   give the successive constraint method
        %   n_field:    number of field variables
        %   n_outputs:  number of outputs
        %   title:      model title
        %   short:      model short description
        %
        % For JRB models the only the model_type and title fields are considered if given, otherwise
        % autodetection is tried.
        Params = [];
        
        % Path to the AffineFunctions file giving the `\theta_i` coefficient functions
        %
        % When using the version 'jrb' (default), leave this field empty in order for the
        % ModelExport to search within the JRB sources for the appropriate AffineFunctions class.
        % There, the package named "models.&lt;modelfolder&gt;" package is searched for the file, where
        % &lt;modelfolder&gt; is the folder name gained by using fileparts on the TargetFolder string.
        %
        % @default ''
        AffFcnsJava = '';

        % Flag that indicates if the export is for a web folder.
        %
        % If true, a file models.txt in the parent folder is searched and opened, and the model
        % directory name added to the list if not already in it. If no such file exists, an attempt
        % is made to create it.
        %
        % If false, no further actions are taken.
        ForWebFolder = true;
    end
    
    methods
        
        function set.TargetFolder(this, value)
            if isempty(value)
                error('Target folder value must not be empty.');
            end
            % Strip trailing slash if set
            if value(end) == '/' || value(end) == '\'
                value = value(1:end-1);
            end
            if exist(value,'dir') ~= 7
                if 1 ~= mkdir(value)
                    error('Could not create directory "%s"', value);
                end
            end
            this.TargetFolder = value;
        end
        
        function set.Version(this, value)
            if ~any(value == [0 1])
                error('Only 0 for "jrb" and 1 for "rbappmit" versions allowed.');
            end
            this.Version = value;
        end
        
        function export_file(this, matfile)
            % function export_file(this, matfile)
            % Exports an m-file containing rb data to a model.
            %
            % Parameters:
            % matfile: The name of the mat-file. Must contain the structs 'model' and
            % 'detailed_data'
            S = load(matfile);
            this.export(S.model, S.detailed_data);
        end
        
        function export(this, model, det)
            % function export(this, model, det)
            % Exports the given model using the detailed data.
            %
            % Parameters:
            % model: The model struct of rbmatlab
            % det: the model's detailed data ('detailed_data')
            
            %% Validity checks
            fprintf('Starting model export to %s ...\n',this.TargetFolder);
            
            %% Re-create the model and reduced data; should not be too costly
            model_data = gen_model_data(model);
            red = gen_reduced_data(model, det);
            
            %% Start the model.xml file with description etc
            f = this.startExportXML(model);
            
            %% Export model data, depending on format
            if this.Version == 0
                this.export_JRB(model, model_data, det, red, f);
            elseif this.Version == 1
                fprintf(f,'<rbappmit_model />\n');
                this.export_rbappmit(model, model_data, det, red);
            end
            
            %% Export affine functions
            this.exportAffineFunctions(f);
            
            %% Export grid (no matter which version - always the rb one for our models)
            this.exportGeometry(model_data.grid, f);
            
            fprintf(f,'\t<visual>\n');
            fprintf(f,'\t\t<plotSteps>%d</plotSteps>\n',model.nt+1);
            fprintf(f,'\t</visual>\n');
            
            this.endExportXML(f);
            
            %% Add folder name to models.txt if it is a web export.
            if this.ForWebFolder
                [parent, foldername] = fileparts(this.TargetFolder);
                mtxt = fullfile(parent, 'models.txt');
                fid = fopen(mtxt,'a+');
                ln = fgetl(fid);
                exists = false;
                while ischar(ln)
                    if strcmp(ln,foldername) 
                        exists = true;
                    end
                    ln = fgetl(fid);
                end
                if ~exists
                    fprintf(fid,'%s\n',foldername);
                end
                fclose(fid);
            end
        end
        
    end 
    methods(Access=private)
        
        function f = startExportXML(this, model)
            % Starts Writing the model.xml file and returns a file handle for further Writing.
            %
            % Already writes the <description> tag.
            %
            f = fopen(fullfile(this.TargetFolder,'model.xml'),'w+');
            fprintf(f,'<?xml version="1.0" encoding="utf-8"?>\n');
            fprintf(f,'<!-- Generated by rbmatlab on %s -->\n',date);
            % Always use big endian format as it's java's native byte order
            version = 'JRB';
            if this.Version == 1
                version = 'rbappmit';
            end
            fprintf(f,'<model type="%s" machformat="be">\n',version);
            fprintf(f,'<description>\n');
            if isfield(this.Params,'title')
                t = this.Params.title;
            else
                t = ['rbmatlab: ' func2str(model.gen_model_data)];
            end
            fprintf(f,'\t<name>%s</name>\n',t);
            fprintf(f,'\t<created>%s</created>\n',datestr(now));
            if isfield(this.Params,'short')
                t = this.Params.short;
            else
                t = ['rbmatlab: ' func2str(model.gen_model_data)];
            end
            fprintf(f,'\t<short>%s</short>\n',t);
            if exist('figure.png','file') == 2
                copyfile('figure.png',fullfile(this.TargetFolder,'figure.png'))
                fprintf(f,'\t<image>figure.png</image>\n');
            end
%             fprintf(f,'\t<infohtml>site_info.html</infohtml>\n');
            fprintf(f,'</description>\n');
        end
        
        function endExportXML(this, f)%#ok
            fprintf(f,'</model>\n');
            fclose(f);
        end
        
        function exportAffineFunctions(this, f)
            % Compile AffFcns class. Classpath needed for JRB model version.
            jfile = this.AffFcnsJava;     
            package = '';
            if this.Version == 0
                if isempty(this.JRBSource)
                    error('For JRB model export version a the JRBSource property must not be empty.');
                end
                if isempty(jfile)
                    [par, fname] = fileparts(this.TargetFolder);
                    package = fullfile('models',fname);
                    jfile = fullfile(fullfile(this.JRBSource,package),'AffineFunctions.java');
                end
                cmd = sprintf('javac -classpath %s -d %s %s',this.JRBSource, this.TargetFolder, jfile);
            elseif this.Version == 1
                if isempty(this.AffFcnsJava)
                    error('Error exporting affine functions: AffFcnsJava field must be set for rbappmit model export.');
                end
                cmd = sprintf('javac -d %s %s',this.TargetFolder, jfile);
            end
            if exist(jfile,'file') ~= 2
                error('Affine functions file "%s" does not exist.',jfile);
            end
            
            fprintf('Compiling and exporting AffineFunctions in "%s"...\n',jfile);
            system(cmd);
            
            outcl = fullfile(fullfile(this.TargetFolder,package),'AffineFunctions.class');
            javajar = fullfile(this.TargetFolder,'classes.jar');
            dexjar = fullfile(this.TargetFolder,'dexclasses.jar');
            
            % Create dex files for Android DalvikVM
            %system(sprintf('dx --dex --output="dexclasses.jar" %s',outcl));
            system(sprintf('dx --dex --no-strict --output="%s" %s',dexjar,outcl));
            
            % Create normal jar file for Java VMs
            system(sprintf('jar -cf %s -C %s %s',javajar,this.TargetFolder,fullfile(package,'AffineFunctions.class')));
            % Remove .class file - we're tidy :-)
            if ~isempty(package)
                rmdir(fullfile(this.TargetFolder,'models'),'s');
                package = ['models.' fname];
            else
                delete(outcl);
            end
            
            if ~isempty(package)
                fprintf(f,'<package>%s</package>\n',package);
            end
        end
        
        function exportGeometry(this, g, f)
            if ~isa(g,'triagrid')
                error('Cannot export an object of class %s as geometry. Instance of triagrid is required.',class(g));
            end
            
            fprintf('Exporting model grid (%s) to geometry files...\n',class(g));
            
            % Write XML entries
            fprintf(f,'<geometry>\n');
            fprintf(f,'\t<dimension>2</dimension>\n');
            fprintf(f,'\t<nodes>%d</nodes>\n',length(g.X));
            fprintf(f,'\t<fieldmapping>ELEMENT</fieldmapping>\n');
            fprintf(f,'</geometry>\n');
            
            % COMMENTED - dont write useless zeros to files, is big enough already.
            % Rather read the model XML and fill up zeros inside the program, if need be.
            % Insert zero z coords (JRB has only 3D vertices)
            %z = single(zeros(size(g.X)));
            
            % Make big column vector
            nodes = single(reshape([g.X g.Y]',[],1));
            
            m = max(g.VI(:));
            if m > intmax('int16')
                error('Too many faces for app: %d. Max allowed is intmax(int16) = %d',m,intmax('int16'));
            end
            faces = int16(reshape(g.VI',[],1));
            
            export.AppExport.saveRealVector(nodes,'vertices.bin',this.TargetFolder);
            export.AppExport.saveRealVector(faces,'faces.bin',this.TargetFolder);
%             save geo nodes faces;
        end
        
        function export_rbappmit(this, model, model_data, det, red)
            %Write the input.in file to define the model problem
            %----------------------------------------------------
            fid = fopen(fullfile(this.TargetFolder,'input.in'),'w');
            
            %check if system type is predefined. If not... make a guess
            %----------------------------------------------------
            if  (isfield(this.Params,'model_type'))
                system_type = this.Params.model_type;
            else
                if isfield(model,'nt')
                    %unsteady...
                    if isfield(model,'affinely_decomposed')&&model.affinely_decomposed
                        system_type = 'LINEAR_UNSTEADY';
                    else
                        system_type = 'NONLINEAR_UNSTEADY';
                    end
                else
                    %steady
                    if isfield(model,'affinely_decomposed')&&model.affinely_decomposed
                        system_type = 'LINEAR_STEADY';
                    else
                        system_type = 'NONLINEAR_STEADY';
                    end
                end
            end
            
            fprintf(fid,'system_type = %s\n',system_type);
            
            %scm type (successive constraint method)
            %---------------------------------------
            if isfield(this.Params,'scm_type')
                scm_type = this.Params.scm_type;
            else
                scm_type = 'NONE';
            end
            fprintf(fid,'scm_type = %s\n',scm_type);
            
            % n_field (number of field variables... bei uns bisher immer nur eine)
            %----------------------------------------------------------------------
            if isfield(this.Params,'n_field')
                n_field = this.Params.n_field;
            else
                n_field = 1; %standard: nur eine Feldvariable
            end
            fprintf(fid,'n_field = %d\n\n',n_field);
            
            %number of components Q
            %----------------------
            fprintf(fid,'Q_a = %d\n',Qa);
            fprintf(fid,'Q_f = %d\n\n',Qf);
            %Kommentar: eigentlich fehlen hier noch Komponenten LL_I des
            %impliziten Operators und die Komponenten der parameter separierten
            %Anfangsbedingungen
            
            %Nmax: maximum number of basis vectors
            %--------------------------------------
            fprintf(fid,'Nmax = %d\n\n',size(det.RB,2));
            
            %Information about parameterization
            %--------------------------------------
            n_parameters = length(model.mu_names);
            fprintf(fid,'n_parameters = %d\n',n_parameters);
            for i=1:n_parameters
                fprintf(fid,['mu',num2str(i-1),'_min = %d\n'],model.mu_ranges{i}(1));
                fprintf(fid,['mu',num2str(i-1),'_max = %d\n'],model.mu_ranges{i}(2));
            end
            
            %Information about outputs
            %-------------------------
            if isfield(this.Params,'n_outputs')
                n_outputs = this.Params.n_outputs;
            else
                n_outputs = 1;
            end
            fprintf(fid,'n_outputs = %d\n',n_outputs);
            
            %Information about time discretization
            %-------------------------------------
            if strcmp(system_type,'LINEAR_UNSTEADY')||strcmp(system_type,'NONLINEAR_UNSTEADY')
                fprintf(fid,'dt = %f\n',(model.T/model.nt)); %Zeitschrittweite
                fprintf(fid,'K = %d\n',model.nt); %number of time steps
                if isfield(model,'theta')
                    fprintf(fid,'euler_theta = %d\n',model.theta);
                else
                    fprintf(fid,'euler_theta = %f\n',0);
                end
            end
            
            %Labels
            %------------------------------------
            if isfield(this.Params,'title')
                fprintf(fid,'title = %s\n',this.Params.title);
            end
            for i=1:n_parameters
                fprintf(fid,['param',num2str(i-1),'_label = %s\n'],model.mu_names{i});
            end
            
            fclose(fid);
            %keyboard
            %-----------------------------------------------------------
            
            %number of basisvectors in RB
            fid = fopen(fullfile(this.TargetFolder,'n_bfs.dat'),'w');
            fprintf(fid,'%d',size(det.RB,2));
            fclose(fid);
            
            %number of nodes in geometry
            fid = fopen(fullfile(this.TargetFolder,'calN.dat'),'w');
            fprintf(fid,'%d',size(model_data.grid.VI,1));
            fclose(fid);
            
            %write output
            disp('Writing output')
            fid=fopen(fullfile(this.TargetFolder,'output_000_000.dat'),'w');
            for i=1:(length(red.s_RB))
                fprintf(fid,'%17.17e ',(red.s_RB(i)));
            end
            fclose(fid);
            
            %write L_E operator components
            disp('Writing LL_E')
            for q=1:Qa
                fid= fopen(fullfile(this.TargetFolder,['RB_A_00',num2str(q-1),'.dat']),'w');
                for i=1:(size(red.LL_E{q},1)^2)
                    fprintf(fid,'%17.17e ',red.LL_E{q}(i));
                end
                fclose(fid);
            end
            
            %write rhs b components
            disp('Writing bb')
            for q=1:Qf
                fid = fopen(fullfile(this.TargetFolder,['RB_F_00',num2str(q-1),'.dat']),'w');
                for i=1:size(red.bb{q},1)
                    fprintf(fid,'%17.17e ', red.bb{q}(i));
                end
                fclose(fid);
            end
            
            %write M matrix (identity matrix)
            disp('Writing M')
            M=eye(size(red.LL_E{1},1), size(red.LL_E{1},2));
            fid = fopen(fullfile(this.TargetFolder,'RB_M_000.dat'),'w');
            for i=1:size(M,1)^2
                fprintf(fid,'%17.17e ',M(i));
            end
            fclose(fid);
            
            %write mass matrices (try=also idendity matrices)
            disp('Writing mass matrices')
            fid = fopen(fullfile(this.TargetFolder,'RB_mass_matrix.dat'),'w');
            for i=1:size(M,1)^2
                fprintf(fid,'%17.17e ',M(i));
            end
            fclose(fid);
            
            fid = fopen(fullfile(this.TargetFolder,'RB_inner_product_matrix.dat'),'w');
            for i=1:size(M,1)^2
                fprintf(fid,'%17.17e ',M(i));
            end
            fclose(fid);
            
            fid = fopen(fullfile(this.TargetFolder,'RB_L2_matrix.dat'),'w');
            for i=1:size(M,1)^2
                fprintf(fid,'%17.17e ',M(i));
            end
            fclose(fid);
            
            %write RB-vectors
            disp('Writing RB-vectors');
            for n=1:size(det.RB,2)
                if n-1<10
                    ns = ['0',num2str(n-1)];
                else
                    ns=num2str(n-1);
                end
                fid = fopen(fullfile(this.TargetFolder,['Z_000_0',ns,'.bin']),'w');
                %for i=1:size(det.RB,1)
                fwrite(fid,det.RB(:,n),'double');
                %end
                fclose(fid);
            end
            
            %write error matrices AqAq
            disp('Writing AqAq')
            Q=length(red.M_EE);
            for q=1:length(red.M_EE)
                fid = fopen(fullfile(this.TargetFolder,['Aq_Aq_00',num2str(floor((q-1)/sqrt(Q))),'_00',num2str(q-floor((q-1)/sqrt(Q))*sqrt(Q)-1),'_norms.bin']),'w');
                fwrite(fid,red.M_EE{q},'double');
                fclose(fid);
            end
            
            %write error matrices FqAq
            disp('Writing FqAq')
            Qb=Qf;
            QL=Qa;
            Fq_Aq=cell(1,QL*Qb);
            for q1=1:Qb
                for q2=1:QL
                    Fq_Aq{q2+(q1-1)*QL} = red.M_Eb{q1+(q2-1)*Qb};
                end
            end
            fid=fopen(fullfile(this.TargetFolder,'Fq_Aq_norms.dat'),'w');
            for i=1:QL*Qb
                for j=1:size(Fq_Aq{i},1)
                    fprintf(fid,'%17.17e ',Fq_Aq{i}(j));
                end
            end
            fclose(fid);
            
            
            %write Aq_Mq norms and Aq_M
            disp('Writing Aq_Mq and Aq_M')
            fid=fopen(fullfile(this.TargetFolder,'Aq_M_norms.dat'),'w');
            for q=1:length(red.M_E)
                for i=1:size(red.M_E{q},1)^2
                    fprintf(fid,'%17.17e ',red.M_E{q}(i));
                end
            end
            fclose(fid);
            fid=fopen(fullfile(this.TargetFolder,'Aq_Mq_norms.dat'),'w');
            for q=1:length(red.M_E)
                for i=1:size(red.M_E{q},1)^2
                    fprintf(fid,'%17.17e ',red.M_E{q}(i));
                end
            end
            fclose(fid);
            
            %write Fq_M and Fq_Mq
            disp('Writing Fq_Mq and Fq_M')
            fid=fopen(fullfile(this.TargetFolder,'Fq_Mq_norms.dat'),'w');
            for q=1:length(red.M_b)
                for i=1:size(red.M_b{q},1)*size(red.M_b{q},2)
                    fprintf(fid,'%17.17e ',red.M_b{q}(i));
                end
            end
            fclose(fid);
            
            fid=fopen(fullfile(this.TargetFolder,'Fq_M_norms.dat'),'w');
            for q=1:length(red.M_b)
                for i=1:size(red.M_b{q},1)*size(red.M_b{q},2)
                    fprintf(fid,'%17.17e ',red.M_b{q}(i));
                end
            end
            fclose(fid);
            
            %write Fq_Fq
            disp('Writing Fq norms');
            fid = fopen(fullfile(this.TargetFolder,'Fq_norms.dat'),'w');
            for q=1:length(red.M_bb)
                fprintf(fid,'%17.17e ',red.M_bb{q});
            end
            fclose(fid);
            
            %write idendity matrices for M_M_norms
            disp('Writing M_M_norms and Mq_Mq_norms');
            fid=fopen(fullfile(this.TargetFolder,'M_M_norms.dat'),'w');
            for i=1:size(M,1)^2
                fprintf(fid,'%17.17e ',M(i));
            end
            fclose(fid);
            
            fid=fopen(fullfile(this.TargetFolder,'Mq_Mq_norms.dat'),'w');
            for i=1:size(M,1)^2
                fprintf(fid,'%17.17e ',M(i));
            end
            fclose(fid);
            
            %output_dual_norm
            disp('Writing dual norm of output')
            fid = fopen(fullfile(this.TargetFolder,'output_000_dual_norms.dat'),'w');
            model.decomp_mode = 0;
            v=model.operators_output(model,model_data);
            fprintf(fid,'%17.17e ',norm(v));
            fclose(fid);
            
            % %write LL_E operator
            % for i=1:Qa
            %     fid = fopen(fullfile(this.TargetFolder,['Aq_Aq_000_00',num2str(i-1),'_norms.bin')],'w');
            %     fwrite(fid,red.LL_E{i},'double');
            %     fclose(fid);
            % end
            %
            %
            % %write mass matrix WRONG!!
            % fid=fopen(fullfile(this.TargetFolder,'RB_inner_product_matrix.dat'),'w');
            % for i=1:(size(model_data.W,1)^2)
            %     fwrite(fid,model_data.W(i),'double');
            %     fwrite(fid,' ');
            % end
            % fclose(fid);
        end
        
        function export_JRB(this, model, model_data, det, red, fxml)
            
            %% Static settings, that might change for different models
            Qm = 1; % number of mass matrices for UNSTEADY-type systems?
            
            %% Determine system type
            if  (isfield(this.Params,'model_type'))
                system_type = this.Params.model_type;
            else
                if isfield(model,'nt')
                    %unsteady...
                    if isfield(model,'affinely_decomposed') && model.affinely_decomposed
                        system_type = 'LINEAR_UNSTEADY';
                    else
                        system_type = 'NONLINEAR_UNSTEADY';
                    end
                else
                    %steady
                    if isfield(model,'affinely_decomposed') && model.affinely_decomposed
                        system_type = 'LINEAR_STEADY';
                    else
                        system_type = 'NONLINEAR_STEADY';
                    end
                end
            end
            
            %% Write rb_model attributes
            % n_field (number of field variables... bei uns bisher immer nur eine)
            %----------------------------------------------------------------------
            if isfield(this.Params,'n_field')
                n_field = this.Params.n_field;
            else
                n_field = 1; %standard: nur eine Feldvariable
            end
            % Number of components Q
            Qa = length(red.LL_E);
            Qf = length(red.bb);
            %Kommentar: eigentlich fehlen hier noch Komponenten LL_I des
            %impliziten Operators und die Komponenten der parameter separierten
            %Anfangsbedingungen
            
            % Nmax: maximum number of basis vectors
            % same as n_bfs, why artificially reduce?
            nmax = size(det.RB,2);
            n_bfs = nmax;       
           
            fprintf(fxml,'<rb_model num_basisfcn="%d" fields="%d" Qa="%d" Qf="%d" Nmax="%d">\n',...
                n_bfs,n_field,Qa,Qf,nmax);
            
            %% System type
            fprintf(fxml,'\t<systype>%s</systype>\n',system_type);
            
            %% Scm type (successive constraint method) --------------------
            if isfield(this.Params,'scm_type')
                scm_type = this.Params.scm_type;
            else
                scm_type = 'NONE';
            end
            fprintf(fxml,'\t<scmtype>%s</scmtype>\n',scm_type);
            
            %% Time discretization ----------------------------------------
            if strcmp(system_type,'LINEAR_UNSTEADY') || strcmp(system_type,'NONLINEAR_UNSTEADY')
                fprintf(fxml,'\t<timeinfo>\n');
                fprintf(fxml,'\t\t<dt>%17.17f</dt>\n',model.T/model.nt);
                fprintf(fxml,'\t\t<K>%d</K>\n',model.nt);
                if isfield(model,'theta')
                    et = model.theta;
                else
                    et = 0;
                end
                fprintf(fxml,'\t\t<euler_theta>%17.17f</euler_theta>\n',et);
                fprintf(fxml,'\t</timeinfo>\n');
            end
            
            %% Information about parameterization -------------------------
            if isfield(model,'mu_names') && ~isempty(model.mu_names)
                pvals = model.get_mu(model);
                n_parameters = length(model.mu_names);
                fprintf(fxml,'\t<parameters number="%d">\n',n_parameters);
                for i=1:n_parameters
                    fprintf(fxml,'\t\t<param name="%s" min="%17.17f" max="%17.17f" label="%s" default="%17.17f"/>\n',...
                        model.mu_names{i},model.mu_ranges{i}(1),model.mu_ranges{i}(2),model.mu_names{i},pvals(i));
                end
                fprintf(fxml,'\t</parameters>\n');
            end
            
            %% Model output
            disp('Writing output');
            export.AppExport.saveRealVector(double(red.s_RB),'output_000_000.bin',this.TargetFolder);
            % Dual norms
            model.decomp_mode = 0;
            v=norm(model.operators_output(model,model_data));
            % The dual norms must equal an Ql*(Ql+1)/2 vector
            % (symmetric information)
            v = this.uppertria(v);
            if numel(v) > 1    
                warning('Model:Export','Export of dual norms for more than one output not yet checked.');
            end
            export.AppExport.saveRealVector(double(v),'output_000_dual_norms.bin',this.TargetFolder);
            
            %% Model initial value data
            disp('Writing initial data');
            for q=1:length(red.a0)
                export.AppExport.saveRealVector(double(red.a0{q}),...
                    sprintf('RB_initial_%.3d.bin',q-1),this.TargetFolder);
            end
            
            
            %% Write reduced basis data
            % L_E operator components (JRB: A_q matrices)
            disp('Writing LL_E');
            for q=1:Qa
                export.AppExport.saveRealMatrix(double(red.LL_E{q}),...
                    sprintf('RB_A_%.3d.bin',q-1),this.TargetFolder);
            end
            
            % Write rhs b components (JRB: F_q vectors)
            disp('Writing bb')
            for q=1:Qf
                export.AppExport.saveRealVector(double(red.bb{q}),...
                    sprintf('RB_F_%.3d.bin',q-1),this.TargetFolder);
            end
            
            % Write M matrix (identity matrix, JRB: M_q matrices)
            disp('Writing M')
            M = double(eye(size(red.LL_E{1},1), size(red.LL_E{1},2)));
            export.AppExport.saveRealMatrix(M, 'RB_M_000.bin',this.TargetFolder);
            
            % JRB: read in TransientRBSystem
            % Here: identity matrix
            export.AppExport.saveRealMatrix(M, 'RB_L2_matrix.bin',this.TargetFolder);
            
            %write RB-vectors
            disp('Writing RB-vectors');
            for n=1:size(det.RB,2)
                export.AppExport.saveRealVector(single(det.RB(:,n)),...
                    sprintf('Z_000_%.3d.bin',n-1),...
                    this.TargetFolder);    
            end
            
            %% Error matrices
            % Write error matrices AqAq
            disp('Writing Aq Aq matrices')
            hlp = this.uppertria(red.M_EE);
            for i=1:Qa
                for j=1:Qa-i+1
                    export.AppExport.saveRealMatrix(double(hlp{(i-1)*Qa + j}),...
                        sprintf('Aq_Aq_%.3d_%.3d_norms.bin',i-1,j-1),...
                        this.TargetFolder);
                end
            end
            
            % Write error matrices Fq Aq (JRB: Fq_Aq representor data in RBSystem)
            % Order looked up in lin_evol_opt_rb_operators
            disp('Writing Fq Aq norms.')
            for i=1:Qa
                for j=1:Qf
                    export.AppExport.saveRealVector(double(red.M_Eb{j+(i-1)*Qf}),...
                        sprintf('Fq_Aq_%.3d_%.3d.bin',j-1,i-1),this.TargetFolder);
                end
            end
            
            %write Fq_Fq
            disp('Writing Fq norms');
            norms = this.uppertria([red.M_bb{:}]);
            export.AppExport.saveRealVector(double(norms),...
                    'Fq_norms.bin',this.TargetFolder); 
            
            % Write Aq_Mq norms
            disp('Writing Aq_Mq')
            for i=1:Qa
                for j=1:Qm
                    export.AppExport.saveRealMatrix(double(red.M_E{(i-1)*Qm+j}),...
                        sprintf('Aq_Mq_%.3d_%.3d_norms.bin',i-1,j-1),this.TargetFolder);
                end
            end
            
            % Write and Fq_Mq
            disp('Writing Fq_Mq')
            for i=1:Qf
                for j=1:Qm
                    export.AppExport.saveRealVector(double(red.M_b{(i-1)*Qm+j}),...
                        sprintf('Fq_Mq_%.3d_%.3d.bin',i-1,j-1),this.TargetFolder);
                end
            end
            
            %write idendity matrices for M_M_norms
            disp('Writing M_M_norms and Mq_Mq_norms');
            for i=1:Qm
                for j=1:Qm-i+1
                    export.AppExport.saveRealMatrix(M,...
                        sprintf('Mq_Mq_%.3d_%.3d.bin',i-1,j-1),this.TargetFolder);
                end
            end
            
            fprintf(fxml,'</rb_model>\n');
        end
        
        function res = uppertria(this, vec)%#ok
            % Returns the upper triangular matrix as linearly addressed vector
            % The input is a vector with a square number of elements,
            % representing a symmetric matrix.
            n = sqrt(numel(vec));
            if n*n ~= numel(vec)
                error('vec has not a square number of elements.');
            end
            res = reshape(vec,n,n);
            res = res(triu(true(n,n)));
        end
    end
    
end