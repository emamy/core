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
    
    properties(SetAccess=private)
                
        DataDirectory;
        
        % The host machine this file model data is created on
        Host;
    end
    
    methods
        function this = FileData(folder)
            if nargin < 1
               folder = KerMor.App.DataStoreDirectory;
            end
            this.DataDirectory = folder;
            this.ensureDir;
            this.Host = KerMor.getHost;
        end
        
        function delete(this)
            fprintf('Called delete on FileData for folder "%s"\n',this.DataDirectory);
            if length(dir(this.DataDirectory)) == 2
                rmdir(this.DataDirectory);
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
                    me = MException('KerMor:data:FileModelData','Could not create dir "%s"',this.DataDirectory);
                    me.addCause(ME);
                    me.throw;
                end
            end
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            % Loads a FileModelData instance.
            %
            % Ensures that the directory associated with this FileData is existent.
            if strcmp(this.host,KerMor.getHost)
                this.ensureDir;
            end
        end
    end
    
end