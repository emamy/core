classdef JaRMoSExport < handle
% JaRMoSExport: Export base class for JaRMoS Models
%
%
% @author Daniel Wirtz @date 2011-08-02
%
% @change{0,6,dw,2012-03-24} Restructured the export functions and classes.
% This class was the former AppExport class.
%
% @change{0,5,dw,2011-10-06} Moved this class to the export package.
%
% @change{0,5,dw,2011-08-25} The saveRealMatrix and saveRealVector now save the matrix in the given
% class, i.e. double,single (=float32, 4bytes), int32 etc.
%
% @new{0,5,dw,2011-08-02} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Constant)
        % The machine format (little/big endian, 32/64bit) to use upon export.
        %
        % @default ieee-be (Big Endian, as Java implements it)
        %
        % See 'doc fopen' for more details.
        MachineFormat = 'ieee-be'; % Changes here also in Line 66!!!
    end
    
     properties
        % The number of degrees of freedom-fields this system computes as output.
        %
        % Real valued and complex valued fields count as one DoF-field.
        %
        % @type integer @default 0
        DoFFields = 0;
        
        % Returns the field descriptions for each logical solution field.
        %
        % Is a struct for each field with the struct-fields
        % - \c Type: The string representation of the jarmos.SolutionFieldType enum
        % - \c Name: The name of the solution field (optional)
        % - \c Mapping: One of the jarmos.geometry.FieldMapping enum string
        % values. (optional, default 'VERTEX')
        %
        % @type struct @default []
        LogicalFields = [];
        
        % A callback for geometry export of the JaRMoS-Model.
        %
        % The function takes two arguments:
        % f: A file handle of the model.xml file ready for writing with fprintf
        % folder: The target folder for binaries
        %
        % @type function_handle @default []
        GeometryExportCallback = [];
        
        % The package of any java classes associated with this model when
        % exported to JaRMoS.
        %
        % @type char @default ''
        JavaExportPackage = '';
        
        % The JaRMoSBase java sources directory. Required if any java classes
        % are to by compiled.
        %
        % @type char @default ''
        JaRMoSBaseSource = '';
        
        % The model type according to the jarmos.ModelType enum
        %
        % @type char @default 'Unknown'
        ModelType = 'Unknown';
        
        % Short description for the model
        %
        % @type char @default ''
        Short = '';
    end
    
    methods
        
        function exportModel(this, model, folder)
            % Exports the given model for use in JaRMoS to the given folder
            %
            % Parameters:
            % model: The model to export. @type models.ReducedModel
            % folder: The target folder @type char
            
            %% Use last folder if KerMor is around
            if nargin == 2
                initialdir = '';
                if ~isempty(which('KerMor'))
                    a = KerMor.App;
                    initialdir = getpref(a.getPrefTag,'LASTDIR','.');
                end
                    folder = uigetdir(initialdir);
                if ~isempty(which('KerMor'))
                    setpref(a.getPrefTag,'LASTDIR',folder);
                end
            end
            
            %% Validations
            if exist(folder,'dir') ~= 7
                try
                    mkdir(folder);
                catch ME
                    mexc = MException('JaRMoSExport:nodir','Directory "%s" did not exist and could not be created.', folder);
                    mexc.addCause(ME);
                    mexc.throw();
                end
            end
            
            %% Create Model XML
            f = fopen(fullfile(folder,'model.xml'),'w+');
            fprintf(f,'<?xml version="1.0" encoding="utf-8"?>\n');
            fprintf(f,'<model type="%s" machformat="be">\n',this.ModelType);
            fprintf(f,'\t<description>\n');
            fprintf(f,'\t\t<short>%s</short>\n',this.Short);
            fprintf(f,'\t\t<created>%s</created>\n',datestr(now,'yyyy-mm-dd'));
            % get image
%             [fn,fo] = rm.createImage;
%             copyfile(fullfile(fo,fn),fullfile(folder,fn));
%             fprintf(f,'\t\t<image>%s</image>\n',fn);
            fprintf(f,'\t</description>\n');
            
            [sourcebase, sourcefiles] = this.typeSpecificExport(f, model, folder);
            
            %% Java sources compilation
            if ~isempty(sourcefiles)
                if isempty(sourcebase)
                    error('No project java sources directory set. Aborting.');
                end
                j = export.Java;
                j.TargetFolder = folder;
                j.Package = this.JavaExportPackage;
                j.JProjectSource = sourcebase;
                j.AdditionalClassPath  = {this.JaRMoSBaseSource};
                j.Sources = sourcefiles;
                j.CreateAndroid = true;
                j.exportFunctions;
            end
            
            fprintf(f,'<numDoFfields>%d</numDoFfields>\n',this.DoFFields);
            
            %% Model geometry
            if ~isempty(this.GeometryExportCallback)
                fprintf(f,'<geometry>\n');
                this.GeometryExportCallback(f, folder, this);
                fprintf(f,'</geometry>\n');
            end
            %% Logical model solution fields
            fi = this.LogicalFields;
            if ~isempty(fi)
                fprintf(f,'<visual><fields>\n');
                for idx=1:length(fi)
                    mapping = 'VERTEX';
                    if isfield(fi(idx),'Mapping') && ~isempty(fi(idx).Mapping)
                        mapping = fi(idx).Mapping;
                    end
                    if ~isempty(fi(idx).Name)
                        fprintf(f,'\t<field type="%s" mapping="%s">%s</field>\n',...
                            fi(idx).Type,mapping,fi(idx).Name);
                    else
                        fprintf(f,'\t<field type="%s" mapping="%s"/>\n',...
                            fi(idx).Type,mapping);
                    end
                end
                fprintf(f,'</fields></visual>\n');
            end
            if ~isempty(this.JavaExportPackage)
                fprintf(f,'<package>%s</package>\n',this.JavaExportPackage);
            end
            fprintf(f,'</model>\n');
            fclose(f);
        end
        
        function saveRealMatrix(this, mat, file, folder)
            % Stores a real double matrix in little endian format.
            %
            % The first 64bits are used to encode the row and column numbers as int32, then the next
            % bytes contain rows*cols*8bytes of double values.
            % The matrix is addressed linearly in a row-wise fashion, i.e. first row, second row ...
            if ~isreal(mat)% || ~ismatrix(mat)
                error('Matrix must contain only real values');
            end
            prec = class(mat);
            % single precision is encoded as float32
            if strcmp(prec,'single')
                prec = 'float32';
            end
            
            if nargin == 3
                if exist(folder,'dir') ~= 7
                    error('Directory "%s" does not exist.', folder);
                end
                file = fullfile(folder,file);
            end
            
            f = fopen(file,'w+',this.MachineFormat);
            try
                [n,m] = size(mat);
                if n > intmax('int32') || m > intmax('int32')
                    error('Cannot save matrix: Dimensions exceed max int32 value.');
                end
                fwrite(f,n,'int32');
                fwrite(f,m,'int32');
%                 for i=1:n
%                     for j=1:m
%                         fwrite(f,mat(i,j),prec);
%                     end
%                 end
                % Matlab Rev. 2009a: entries are written in column order, which must be transposed
                % in order to read them in row order (native storage format in java
                % apache.commons.math
                fwrite(f,mat',prec);
            catch ME
                fclose(f);
                rethrow(ME);
            end
            fclose(f);
        end
        
        function saveRealVector(this, vec, file, folder)
            % Stores a real double vector in little endian format.
            %
            % The first 64bits are used to encode the vector size as int32, then the next
            % bytes contain size*8bytes of double values.
            if ~isreal(vec) || ~isvector(vec)
                error('vec must be a vector with real values');
            end
            prec = class(vec);
            % single precision is encoded as float32
            if strcmp(prec,'single')
                prec = 'float32';
            end
            
            if nargin == 3
                if exist(folder,'dir') ~= 7
                    error('Directory "%s" does not exist.', folder);
                end
                file = fullfile(folder,file);
            end
            
            f = fopen(file,'w+',this.MachineFormat);
            try
                s = length(vec);
                if s > intmax('int32')
                    error('Cannot save vector: Dimension exceeds max int32 value.');
                end
                fwrite(f,s,'int32');
%                 for i=1:s
%                     fwrite(f,vec(i),prec);
%                 end
                fwrite(f, vec, prec);
            catch ME
                fclose(f);
                rethrow(ME);
            end
            fclose(f);
        end
    end
    
    methods(Abstract, Access=protected)
        [sourcebase, sourcefiles] = typeSpecificExport(this, f, model, settings, folder);
    end
    
end