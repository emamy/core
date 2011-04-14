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
    % @new{0,3,dw,2011-04-12} Added this class to aid new class creation and offer customized class
    % / function skeletons.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Dependent)
        % The current developing author
        Author;
        
        % Author's short tag for identification in doxygen @@new and @@change tags
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
            inst = Devel;
        end
        
        function setup
            d = Devel.Instance;
            fprintf('Running KerMor developer setup...\n');
            str = sprintf('Please enter your full name: ');
            d.Author = input(str,'s');
            str = sprintf('Please enter your author''s short tag: ');
            d.AuthorShort = input(str,'s');
            fprintf('Finished!\n');
        end
        
        function newClass(name, overwrite)
            % Creates a new class using the template_class file
            %
            % Parameters:
            % name: The target class name
            % overwrite: [Optional] Overwrite any existing file
            %
            % @todo check class name validity!
            if ~isa(name,'char')
                error('The class name must be a string');
            end
            if nargin < 2
                overwrite = false;
            end
            Devel.process(name, 'template_class', overwrite);
        end
        
        function newFun(name, overwrite)
            % Creates a new function using the template_fun file
            %
            % Parameters:
            % name: The target function name
            % overwrite: [Optional] Overwrite any existing file
            %
            % @todo check fun name validity!
            if nargin < 2
                overwrite = false;
            end
            
            if ~isa(name,'char')
                error('The class name must be a string');
            end
            
            Devel.process(name, 'template_fun', overwrite);
        end
        
    end
    
    methods(Static,Access=private)
        function process(name, template_file, overwrite)
            dir = pwd;
            fname = fullfile(dir,[name '.m']);
            if ~overwrite && exist(fname,'file')
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

