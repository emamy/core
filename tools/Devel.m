classdef Devel < handle
    % Developer utilities.
    %
    % New function/class generation:
    % In order to create new classes/functions with a customized initial comment etc. use the static
    % methods Devel.newClass and Devel.newFun functions. They use the template files template_class
    % and template_fun from the KerMor root directory and produce the new class by replacing any of
    % the following occurences by the respective same-named properties of any Devel instance:
    % - \c $author The full name of the author
    % - \c $authorshort The author's short tag for identification in the documentation
    % - \c $mainver The current KerMor main version
    % - \c $subver The current KerMor sub version
    % - \c $date Today's date in the form yyyy-mm-dd
    %
    % See the template files for examples. 
    %
    % @author Daniel Wirtz @date 2011-04-12
    %
    % @change{0,7,dw,2013-01-23} Removed the overwrite flag and added option to specify a path
    % in which the file should be created.
    %
    % @change{0,4,sa,2011-06-02} Added comments to the function setup
    %
    % @new{0,3,dw,2011-04-12} Added this class to aid new class creation and offer customized class
    % / function skeletons.
    %
    % @todo Check class/function name string validity for newClass/newFun
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(Dependent)
        % The current developing author
        %
        % @type char
        Author;
        
        % Author's short tag for identification in doxygen @@new and @@change tags
        % @type char
        AuthorShort;
        
        % Always today's date
        Date;
    end
    
    methods
        function set.Author(this, value)%#ok
            setpref('KERMOR_DEVEL','author',value);
        end
        function set.AuthorShort(this, value)%#ok
            setpref('KERMOR_DEVEL','authorshort',value);
        end
        
        function v = get.Author(this)%#okdate
            v = getpref('KERMOR_DEVEL','author','<No author set>');
        end
        function v = get.AuthorShort(this)%#ok
            v = getpref('KERMOR_DEVEL','authorshort','<No author short set>');
        end
        function v = get.Date(this)%#ok
            v = datestr(now, 'yyyy-mm-dd');
        end
    end
    
    methods(Static)
        
        function inst = Instance
            % Returns the singleton instance of Devel.
            %
            % Return values:
            % inst: The instance @type Devel
            inst = Devel;
        end 
        
        function setup
            % Setup variables for Kermor Development when Kermor is setup
            % for the first time or any other time manually
            %
            % Add the name of the author and his initials(short tag)
            d = Devel.Instance;
            fprintf('Running KerMor developer setup...\n');
            str = sprintf('Please enter your full name: ');
            d.Author = input(str,'s');
            str = sprintf('Please enter your author''s short tag: ');
            d.AuthorShort = input(str,'s');
            fprintf('Finished!\n');
        end
        
        function newClass(name, varargin)
            % Creates a new class using the template_class file
            %
            % Parameters:
            % name: The target class name
            % varargin: Additional parameters.
            % dir: The target directory @type char @default pwd
            if ~isa(name,'char')
                error('The class name must be a string');
            end
            Devel.process(name, 'template_class', varargin{:});
        end
        
        function newFun(name, varargin)
            % Creates a new function using the template_fun file
            %
            % Parameters:
            % name: The target function name
            % varargin: Additional parameters.
            % dir: The target directory @type char @default pwd
            %
            % @todo check fun name validity!
            
            if ~isa(name,'char')
                error('The class name must be a string');
            end
            
            Devel.process(name, 'template_fun', varargin{:});
        end
        
    end
    
    methods(Static,Access=private)
        function process(name, template_file, dir)
            if nargin < 3
                dir = pwd;
            end
            fname = fullfile(dir,[name '.m']);
            if exist(fname,'file')
                error('File %s already exists.',fname);
            end
            
            % Open template
            fh = fopen(template_file);
            str = fscanf(fh,'%c');
            fclose(fh);
            
            % Replace occurences
            %str = strrep(str,'%','%%')
            d = Devel.Instance;
            str = strrep(str,'$name',name);
            str = strrep(str,'$date',d.Date);
            str = strrep(str,'$authorshort',d.AuthorShort);
            str = strrep(str,'$author',d.Author);
            str = strrep(str,'$mainver',KerMor.MainVersion);
            str = strrep(str,'$subver',KerMor.SubVersion);
            
            % Write to file
            fh = fopen(fname, 'w+');
            fprintf(fh, '%c', str);
            fclose(fh);
            
            % Open for editing
            edit(fname);
        end
    end
    
end

