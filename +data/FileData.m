classdef FileData < handle
% FileData: 
%
%
%
% @author Daniel Wirtz @date 2012-07-09
%
% @new{0,6,dw,2012-07-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=protected)
        % The root folder where the FileData's files are saved.
        %
        % @type char @default KerMor.DataStoreDirectory
        DataDirectory;
        
        % The host machine this file data is created on.
        %
        % If the 
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
            if nargin < 1
               error('A root folder must be specified.');
            end
            this.DataDirectory = data_dir;
            this.ensureDir;
            this.Host = KerMor.getHost;
        end
        
        function relocate(this, newDataDirectoryRoot)%#ok
            % @TODO Implement multiple-host FileData's with different storage_root folders (per
            % host for example)
        end
        
        function delete(this)
            if ~this.isSaved && exist(this.DataDirectory,'dir') == 7
                if length(dir(this.DataDirectory)) ~= 2
                    warning('KerMor:FileData',...
                        'A FileData instance (%s) should be deleted but the DataDirectory "%s" is not empty. Not deleting.',...
                        class(this),this.DataDirectory);
                else
                    rmdir(this.DataDirectory);
                end
            end
        end
    end
    
    methods(Access=protected)
        function file = getfile(this, file)
            file = fullfile(this.DataDirectory,file);
            if exist(file,'file') ~= 2
                error('File "%s" not found in folder "%s". Have you deleted data files?',file,this.DataDirectory);
            end
        end
        
        function ensureDir(this)
            if exist(this.DataDirectory,'dir') ~= 7
                try
                    mkdir(this.DataDirectory);
                catch ME
                    me = MException('KerMor:FileData','Could not create directory "%s"',this.DataDirectory);
                    me.addCause(ME);
                    me.throw;
                end
            end
        end
    end
    
    methods(Access=protected)
        function this = saveobj(this)
            % Set saved flag so that the data files do not get deleted on the delete method
            this.isSaved = true;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            % Loads a FileData instance.
            %
            % Ensures that the directory associated with this FileData is existent.
            if strcmp(this.Host,KerMor.getHost)
                this.ensureDir;
            end
        end
    end
    
end