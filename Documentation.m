classdef Documentation
    % Documentation: Class for information & tasks concerning the KerMor
    % documentation.
    %
    %
    % @author Daniel Wirtz @date 2011-10-13
    %
    % @new{0,5,dw,2011-10-13} Added this class and moved documentation related
    % stuff here from the KerMor class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Constant)
        % The KerMor documentation directory, i.e. where createDocs places
        % the generated documentation.
        %
        % Can be set during Documentation.setup
        % Readonly.
        DocumentationDirectory = 'df';%getenv('KERMOR_DOCS');
        
        % The doxygen binary used to create the documentation.
        %
        % Can be set during Documentation.setup
        % Readonly.
        Doxygen = 'asdf';%getenv('KERMOR_DOXYBIN');
    end
    
    properties(Dependent)
        % Returns where the documentation is located.
        %
        % This is either a web-site or the local documentation.
        DocumentationLocation;
    end
    
    methods
        function d = get.DocumentationLocation(this)%#ok
            d = getenv('KERMOR_DOCS');
            if isempty(d) || ~exist(fullfile(d,'index.html'),'file')
                d = 'http://www.agh.ians.uni-stuttgart.de/documentation/kermor';
            end
        end
    end
    
    methods(Static)
        
        function create(uml, open)
            % Creates the Doxygen documentation
            %
            % Parameters:
            % uml: Set to true to create UML-like graphics output @default false @type boolean
            % open: Set to true if the documentation should be opened after
            % successful compilation @default false @type bool
            if nargin < 2
                open = false;
                if nargin < 1
                    uml = false;
                end
            end
            
            %% Operation-system dependent actions
            if isunix
                Documentation.createUnix(uml, open);
            elseif ispc
                Documentation.createWindows(uml, open);
            end
        end
        
        function setup
            % Performs the setup of environment variables needed to create
            % the KerMor documentation.
            %
            % Can be called by its own, however, typically this method is
            % called automatically during the KerMor.install script.
            fid = fopen('~/.bashrc','a+');
            try
                %% Documentation directory
                if isempty(getenv('KERMOR_DOCS'))
                    d = fullfile(h,'documentation','output');
                    str = sprintf(['No documentation output directory has been set yet.\n'...
                        'The default will be %s\n'...
                        'Do you want to specify a custom output directory? (Y)es/(N)o: '],d);
                    ds = lower(input(str,'s'));
                    if isequal(ds,'y')
                        d = uigetdir(h,'Please select the documentation output folder.');
                        if d == 0
                            d = fullfile(h,'documentation','output');
                            fprintf('Operation cancelled, using default directory %s...\n',d);
                        end
                    end
                    fprintf(fid,'export KERMOR_DOCS="%s"\n',d);
                    setenv('KERMOR_DOCS',d)
                end
                
                %% Doxygen binary
                if isempty(getenv('KERMOR_DOXYBIN'))
                    [s,r] = system('which doxygen');
                    db = [];
                    if ~isempty(r)
                        db = 'doxygen';
                        str = sprintf(['Doxygen installation is available (%s).\n'...
                            'Do you want to use a different doxygen binary? (Y)es/(N)o: '],strrep(r,char(10),''));
                        yn = lower(input(str,'s'));
                    end
                    if isempty(r) || isequal(yn,'y')
                        [f,p] = uigetfile('~/*.*','Select the custom doxygen binary file.');
                        if f ~= 0
                            db = fullfile(p, f);
                        end
                    end
                    if ~isempty(db)
                        fprintf(fid,'export KERMOR_DOXYBIN="%s"\n',db);
                        setenv('KERMOR_DOXYBIN',db);
                    else
                        warning('Documentation:setup','No doxygen binary selected. Documentation creation will not work.');
                    end
                end
                
                fclose(fid);
            catch ME
                fclose(fid);
                rethrow(ME);
            end
        end
    end
    
    methods(Static,Access=private)
        function createUnix(uml, open)
            % Creates the KerMor documentation on UNIX paltforms.
            %
            
            cmd = fullfile(KerMor.App.HomeDirectory,'documentation','make.sh');
            % Add version argument
            cmd = [cmd ' ' KerMor.MainVersion '.' KerMor.SubVersion];
            % Add uml argument
            if uml
                cmd = [cmd ' uml'];
            end
            [s,r] = system(cmd);
            wpos = strfind(r,'Logged warnings:');
            if ~isempty(wpos)
                endpos = strfind(r,'Complete log file');
                cprintf([0 .5 0],r(1:wpos-1));
%                 cprintf([1,.4,0],strrep(r(wpos:endpos-1),'\','\\'));
                cprintf([0 .5 0],r(endpos:end));
            else
%                 cprintf([0 .5 0],strrep(r,'\','\\'));
            end
            fprintf('\n');
            index = fullfile(getenv('KERMOR_DOCS'), 'index.html');
            if open
                % Try to use iceweasel per default (nasty, uhm?)
                [s,r] = system('which iceweasel');
                if ~isempty(r)
                    cmd = 'iceweasel ';
                else
                    % Otherwise: use user preferred browser
                    cmd = 'xdg-open ';
                end
                system([cmd index]);
            end
        end
        
        function createWindows(uml, open)%#ok
            % Creates the documentation on a windows platform.
            %
            % Currently not implemented.
            error('Creating documentation on Windows is not yet implemented.');
        end
    end
    
end