classdef FileData < handle
% FileData: Base class for access of files stored in a specific folder in the local file system
%
% @author Daniel Wirtz @date 2012-07-09
%
% @new{0,6,dw,2012-07-09} Added this class.
%
% @change{0,7,dw,2013-03-27} Improved the modularity of FileData with more verbosity and
% implemented a basic "relocate" method.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=protected)
        % The root folder where the FileData's files are saved.
        %
        % @type char @default ''
        DataDirectory = '';
        
        % The host machine this file data is created on.
        %
        % @type char @default ''
        Host = '';
    end
    
    properties(Access=protected)
        % This flag indicates that this FileData instance has been stored to disk via the save
        % method somewhere.
        %
        % In subclasses, this can be used to clean up any files in the DataDirectory when the
        % object's delete method is called.
        %
        % @type logical @default false
        %
        % See also: DataDirectory delete
        isSaved = false;
    end
    
    methods
        function this = FileData(data_dir)
            if nargin == 1
                this.DataDirectory = data_dir;
            end
            this.Host = KerMor.getHost;
        end
        
        function relocate(this, newDataDirectoryRoot)
            % Relocates this FileData instance to a different folder.
            %
            % Currently, moving it's files must be done manually.
            %
            % Parameters:
            % newDataDirectoryRoot: The new data directory. @type char
            %
            % @todo Implement multiple-host FileData's with different storage_root folders (per
            % host for example)
            this.DataDirectory = newDataDirectoryRoot;
            this.Host = KerMor.getHost;
        end
        
        function delete(this)
            if ~isempty(this.DataDirectory) && ~this.isSaved && exist(this.DataDirectory,'dir') == 7
                if length(dir(this.DataDirectory)) ~= 2
                    warning('KerMor:FileData',...
                        'A FileData instance (%s) should be deleted but the DataDirectory "%s" is not empty. Not deleting.',...
                        class(this),this.DataDirectory);
                else
                    rmdir(this.DataDirectory);
                end
            end
        end        
        
        function set.DataDirectory(this, value)
            if ~isempty(value)
                if Utils.ensureDir(value);
                    this.DataDirectory = value;
                    return;
                else
                    warning('KerMor:FileData','Could not make sure that the directory "%s" exists.\nPlease ensure a correct path and try again.',value);
                end    
            end
            this.DataDirectory = '';    
        end
    end
    
    methods(Access=protected)
        function file = getfile(this, file)
            file = fullfile(this.DataDirectory,file);
            if exist(file,'file') ~= 2
                if strcmp(this.Host,KerMor.getHost)
                    error('File "%s" not found in folder "%s". Have you deleted data files?',...
                        file,this.DataDirectory);
                else
                    error(['File "%s" not found. This FileData instance was created on host "%s"'...
                        'and the DataDirectory cannot be found on the local machine "%s".\n'...
                        'Copy the files to an accessible directory and use the "relocate" method.'],...
                        file,this.Host,KerMor.getHost);
                end
            end
        end
        
        function this = saveobj(this)
            % Set saved flag so that the data files do not get deleted on the delete method
            this.isSaved = true;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this, initfrom)
            % Loads a FileData instance.
            %
            % Ensures that the directory associated with this FileData is existent.
            
            created = false;
            if ~isa(this, 'data.FileData')
                initfrom = this;
                this = data.FileData(initfrom.DataDirectory);
                created = true;
            end
            if nargin == 2 || created
                this.DataDirectory = initfrom.DataDirectory;
                this.Host = initfrom.Host;
                this.isSaved = initfrom.isSaved;
            end
            if ~isempty(this.DataDirectory) && ~Utils.ensureDir(this.DataDirectory) 
                str = 'Could not create/access "%s". Loading files will not work until fixed.';
                if ~strcmp(this.Host,KerMor.getHost)
                    str = sprintf(['%s\nNote that this FileData has been created on host %s.'...
                        'If you changed the host (local: %s), use the "relocate" method.'],str,this.Host,KerMor.getHost);
                end
                warning('KerMor:FileData',str,this.DataDirectory);
            end
        end
    end
    
end