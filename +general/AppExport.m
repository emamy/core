classdef AppExport
% AppExport: Export class for Android KerMORDSApp
%
%
% @author Daniel Wirtz @date 2011-08-02
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
        MachineFormat = 'ieee-be';
    end
    
    methods(Static)
        
        function exportReducedModel(rm, folder)
            %% Validations
            if exist(folder,'dir') ~= 7
                try
                    mkdir(folder);
                catch ME
                    m = MException('KerMor:AppExport:nodir','Directory "%s" did not exist and could not be created.', folder);
                    m.addCause(ME);
                    m.throw();
                end
            end
            
            %% Create Model XML
            f = fopen(fullfile(folder,'model.xml'),'w+');
            fprintf(f,'<?xml version="1.0" encoding="utf-8"?>\n');
            fprintf(f,'<model type="kermor" machformat="be">\n');
            fprintf(f,'\t<description>\n');
            fprintf(f,'\t\t<short>%s</short>\n',rm.Name);
            % get image
            [fn,fo] = rm.createImage;
            copyfile(fullfile(fo,fn),fullfile(folder,fn));
            fprintf(f,'\t\t<image>%s</image>\n',fn);
            
            fprintf(f,'\t</description>\n');
            fprintf(f,'\t<data><kermor>\n');
            fprintf(f,'\t\t<T>%17.17f</T>\n',rm.T);
            fprintf(f,'\t\t<dt>%17.17f</dt>\n',rm.dt);
            
            % Get system dimension from initial value
            mu = [];
            if ~isempty(rm.ParamSamples)
                mu = rm.ParamSamples(:,1);
            end
            fprintf(f,'\t<dim>%d</dim>\n',size(rm.System.x0.evaluate(mu),1));
            
            %% Export model data
            general.AppExport.saveRealMatrix(rm.V,'V.bin',folder);
            general.AppExport.saveRealMatrix(rm.W,'W.bin',folder);
            %general.AppExport.saveRealMatrix(rm.W,'G.bin',folder);
            
            %% Export system data
            
            % Parameters
            fprintf(f,'\t<parameters>\n');
            p = rm.FullModel.System.Params;
            if ~isempty(p)
                %pvals = zeros(length(p),2);
                for k=1:length(p)
                    fprintf(f,'\t\t<param name="%s" min="%17.17f" max="%17.17f" label="%s"/>\n',...
                        p(k).Name,p(k).MinVal,p(k).MaxVal,p(k).Name);
                end
                %general.AppExport.saveRealMatrix(pvals, 'paramvalues.bin', folder);
            end
            fprintf(f,'\t</parameters>\n');
            
            % Kernel expansion
            s = rm.System;
            cf = s.f;
            fprintf(f,'\t<corefun>\n');
            fprintf(f,'\t\t<corefuntype>%s</corefuntype>\n',class(cf));
            if isa(cf,'kernels.KernelExpansion')
                general.AppExport.saveRealMatrix(cf.Ma,'Ma.bin',folder);
                general.AppExport.saveRealMatrix(cf.Centers.xi,'xi.bin',folder);
                exportKernel(cf.Kernel,'kernel.bin',folder);
                fprintf(f,'\t\t<statekernel>%s</statekernel>\n',class(cf.Kernel));
                if isa(cf,'kernels.ParamTimeKernelExpansion')
                    general.AppExport.saveRealVector(cf.Centers.ti,'ti.bin',folder);
                    exportKernel(cf.TimeKernel,'timekernel.bin',folder);
                    fprintf(f,'\t\t<timekernel>%s</timekernel>\n',class(cf.TimeKernel));
                    
                    general.AppExport.saveRealMatrix(cf.Centers.mui,'mui.bin',folder);
                    exportKernel(cf.ParamKernel,'paramkernel.bin',folder);
                    fprintf(f,'\t\t<paramkernel>%s</paramkernel>\n',class(cf.ParamKernel));
                end
            else
                error('System function must be a kernel expansion.');
            end
            fprintf(f,'\t</corefun>\n');
            
            % Input
            if ~isempty(s.B)
                if isa(s.B,'dscomponents.LinearInputConv')
                    general.AppExport.saveRealMatrix(s.B.B,'B.bin',folder);
                    
                elseif isa(s.B,'dscomponents.AffLinInputConv')
                    
                end
                fprintf(f,'\t<inputconvtype>%s</inputconvtype>\n',class(s.B));
            end
            
            % Output
            if ~isempty(s.C)
                if isa(s.C,'dscomponents.LinearOutputConv')
                    general.AppExport.saveRealMatrix(s.C.C,'C.bin',folder);
                elseif isa(s.B,'dscomponents.AffLinOutputConv')
                    
                end
                fprintf(f,'\t<outputconvtype>%s</outputconvtype>\n',class(s.C));
            end
            
            % Initial value
            if isa(s.x0,'dscomponents.ConstInitialValue')
                general.AppExport.saveRealVector(s.x0.x0,'x0.bin',folder);
            elseif isa(s.B,'dscomponents.AffineInitialValue')
                
            end
            fprintf(f,'\t<initialvaluetype>%s</initialvaluetype>\n',class(s.x0));
            
            fprintf(f,'\t</kermor></data>\n');
            fprintf(f,'</model>\n');
            fclose(f);
            
            function exportKernel(k, file, folder)
                if isa(k,'kernels.GaussKernel')
                    general.AppExport.saveRealVector(k.Gamma,file,folder);
                elseif isa(k,'kernels.LinearKernel')
                    % do nothing.
                end
            end
            
        end
        
        function saveRealMatrix(mat, file, folder)
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
            
            f = fopen(file,'w+',general.AppExport.MachineFormat);
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
        
        function saveRealVector(vec, file, folder)
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
            
            f = fopen(file,'w+',general.AppExport.MachineFormat);
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
    
end