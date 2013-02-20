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
        % The names of the source files.
        %
        % If a JProjectSource is given, just specify the filenames. It will
        % be used in combination with a possibly set package to determine
        % the source file location.
        %
        % @type cell<char>
        Sources;
        
        % An additional JProjectSource for inclusion.
        %
        % Set this value if the class implements Interfaces from their
        % respective Java projects.
        %
        % @type char @default ''
        JProjectSource = '';
        
        % An additional cell of paths for the java classpath.
        %
        % Set this value if the class implements Interfaces from their
        % respective Java projects.
        %
        % @type cell<char> @default {}
        AdditionalClassPath = {};
        
        % The classes package (optional, required if it has one but the default package)
        %
        % @type char @default ''
        Package = '';
        
        % The target folder for the class output.
        %
        % Required, an error will be thrown if the folder cannot be
        % created.
        %
        % @type char @default ''
        TargetFolder = '';
        
        % Create a class for java virtual machine
        %
        % @default true @type logical
        CreateJVM = true;
        
        % Create a class for the Dalvik VM (Android dex required)
        %
        % @default true @type logical
        CreateAndroid = false;
    end
    
    properties(Constant)
        JavaJarFile = 'classes.jar';
        
        AndroidJarFile = 'dexclasses.jar';
    end
    
    methods
        
        function exportFunctions(this)
            % Compile AffFcns class. JProjectSource needed for JRB/JKerMor models.
            if isempty(this.Sources)
                error('Sourcefile must be set');
            end
            
            outjava = cell.empty(1,0);
            outdex = cell.empty(1,0);
            for sidx = 1:length(this.Sources)
                src = this.Sources{sidx};
                pkgdir = strrep(this.Package,'.',filesep);
                if isempty(this.JProjectSource)
                    cmd = sprintf('javac -d %s %s', this.TargetFolder, [src '.java']);
                else
                    jfile = fullfile(fullfile(this.JProjectSource,pkgdir),[src '.java']);
                    cp = this.JProjectSource;
                    if ~isempty(this.AdditionalClassPath)
                        cp = Utils.implode([{cp}, this.AdditionalClassPath],pathsep,'%s');
                    end
                    cmd = sprintf('javac -classpath "%s" -d %s %s', cp, this.TargetFolder, jfile);
                end

                fprintf('Compiling and exporting class "%s"...\n',src);
                system(cmd);
                
                relpath = fullfile(pkgdir,[src '.class']);
                if this.CreateJVM
                    outjava{sidx} = sprintf('-C %s %s',this.TargetFolder,relpath);
                end
                if this.CreateAndroid
                    outdex{sidx} = fullfile(this.TargetFolder,relpath);
                end
            end
            
            % Create normal jar file for Java VMs
            if this.CreateJVM
                javajar = fullfile(this.TargetFolder,this.JavaJarFile);
                system(sprintf('jar -cf %s %s',javajar,...
                    Utils.implode(outjava,' ','%s')));
            end
            % Create dex files for Android DalvikVM
            if this.CreateAndroid
                dexjar = fullfile(this.TargetFolder,this.AndroidJarFile);
                %system(sprintf('dx --dex --output="dexclasses.jar" %s',outcl));
                system(sprintf('dx --dex --no-strict --output="%s" %s',dexjar,...
                    Utils.implode(outdex,' ','%s')));
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
        
        function set.TargetFolder(this, value)
            if (isunix && value(1) ~= '/') || ispc && isempty(strfind(value,':'))
                value = fullfile(pwd,value);
            end
            this.TargetFolder = value;
        end
        
    end
    
end