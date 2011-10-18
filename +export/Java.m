classdef Java < handle
% Java: Java utils like compiling classes out of matlab.
%
%
%
% @author Daniel Wirtz @date 2011-09-23
%
% @new{0,5,dw,2011-09-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The name of the source file.
        %
        % If a JProjectSource is given, just specify the filename. It will
        % be used in combination with a possibly set package to determine
        % the source file location. 
        Sourcefile;
        
        % An additional JProjectSource for inclusion.
        %
        % Set this value if the class implements Interfaces from their
        % respective Java projects.
        %
        % @default ''
        JProjectSource = '';
        
        % The classes package (optional, required if it has one but the default package)
        Package = '';
        
        % The target folder for the class output.
        %
        % Required, an error will be thrown if the folder cannot be
        % created.
        %
        % @default ''
        TargetFolder;
        
        % Create a class for java virtual machine
        %
        % @default true
        CreateJVM = true;
        
        % Create a class for the Dalvik VM (Android dex required)
        %
        % @default true
        CreateAndroid = false;
    end
    
    properties(Constant)
        JavaJarFile = 'classes.jar';
        
        AndroidJarFile = 'dexclasses.jar';
    end
    
    methods
        
        function exportFunctions(this)
            % Compile AffFcns class. JProjectSource needed for JRB/JKerMor models.
            if isempty(this.Sourcefile)
                error('Sourcefile must be set');
            end
            
            [dummy, fname] = fileparts(this.Sourcefile);
            pkgdir = strrep(this.Package,'.',filesep);
            if isempty(this.JProjectSource)
                cmd = sprintf('javac -d %s %s', this.TargetFolder, this.Sourcefile);
            else
                jfile = fullfile(fullfile(this.JProjectSource,pkgdir),this.Sourcefile);
                cmd = sprintf('javac -classpath %s -d %s %s', this.JProjectSource, this.TargetFolder, jfile);
            end

            fprintf('Compiling and exporting class "%s"...\n',this.Sourcefile);
            system(cmd);
            
            outcl = fullfile(fullfile(this.TargetFolder,pkgdir),[fname '.class']);
            
            % Create normal jar file for Java VMs
            if this.CreateJVM
                javajar = fullfile(this.TargetFolder,this.JavaJarFile);
                system(sprintf('jar -cf %s -C %s %s',javajar,this.TargetFolder,fullfile(pkgdir,[fname '.class'])));
            end
            % Create dex files for Android DalvikVM
            if this.CreateAndroid
                dexjar = fullfile(this.TargetFolder,this.AndroidJarFile);
                %system(sprintf('dx --dex --output="dexclasses.jar" %s',outcl));
                system(sprintf('dx --dex --no-strict --output="%s" %s',dexjar,outcl));
            end
            
            % Remove .class file - we're tidy :-)
            if ~isempty(this.Package)
                folder = this.Package(1:strfind(this.Package,'.')-1);
                if isempty(folder)
                    folder = this.Package;
                end
                rmdir(fullfile(this.TargetFolder,folder),'s');
            else
                delete(outcl);
            end           
        end
        
%         function set.Sourcefile(this, value)
%             if exist(value,'file') ~= 2
%                 error('Source file "%s" does not exist.',value);
%             end
%             this.Sourcefile = value;
%         end

        function set.CreateJVM(this, value)
            if system('javac') == 1
                error('Error checking for javac compiler.');
            end
            this.CreateJVM = value;
        end
        
        function set.CreateAndroid(this, value)
%             if system('dx') == 1
%                 error('Error checking for dex compiler.');
%             end
            this.CreateAndroid = value;
        end
        
    end
    
end