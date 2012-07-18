classdef FileTrajectoryData < data.ATrajectoryData & data.FileData & general.ABlockSVD
% FileTrajectoryData: Trajectory data stored in external files.
%
% The constructor takes an optional storage_root parameter.
% If given, it must be a valid folder in the file system.
% If not given, the value set in KerMor.App.DataStoreDirectory is used.
%
% @author Daniel Wirtz @date 2011-08-04
%
% @change{0,6,dw,2012-04-27} Added an internal property 'host' to the class
% in order to record on which machine the FileTrajectoryData instance was
% created on. When saved and loaded at another machine, the dictionary will
% be cleared if the same path does not exist anymore on the new machine.
%
% @new{0,5,dw,2011-11-02} Implemented the getBoundingBox method from superclass. Now keeping
% track of bounding box while adding trajectories.
%
% @new{0,5,dw,2011-10-14} Added a new consolidate method in order to
% rebuild the internal index hashmap from the files in a directory and the
% data.ATrajectoryData.ParamSamples. This method is called after parallel
% trajectory computation as the remotely constructed hashmaps are not joint
% back together as one.
%
% @new{0,5,dw,2011-08-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        % The HashMap used to store the indices for each trajectory.
        hm;
        
        % Stores the DoFs of the trajectories and parameter sizes
        sizes = [];
        
        % Bounding box minimum
        bbmin = [];
        
        % Bounding box maximum
        bbmax = [];
        
        trajlen = [];
    end
    
    methods
        function this = FileTrajectoryData(varargin)
            % Creates a new ModelData instance with trajectory DataDirectory in a file folder.
            %
            % Parameters:
            % varargin: Either a data.ModelData instance to infer the storage root from, or a
            % string containing a valid folder. @default A temporary folder within the
            % KerMor.TempDirectory
            if ~usejava('jvm')
                error('FileTrajectoryData cannot be used as java is not enabled.');
            end
            if isempty(varargin)
                data_dir = fullfile(KerMor.App.TempDirectory,sprintf('temp_ftd_%s',...
                    general.IDGenerator.generateID));
            elseif isa(varargin{1},'data.ModelData')
                data_dir = fullfile(varargin{1}.DataDirectory,'trajectories');
            elseif ischar(varargin{1})
                data_dir = varargin{1};
            else
                error('Invalid argument: %s',class(varargin{1}));
            end
            this = this@data.FileData(data_dir);
            this.hm = java.util.HashMap;
            this.clearTrajectories;
        end
        
        function delete(this)
            % Destructor for FileTrajectoryData
            %
            % Deletes the DataDirectory if no trajectories are stored in it or it has not been
            % saved somewhere.
            if ~this.isSaved || this.hm.size == 0
                this.clearTrajectories;
            end
            % Superclass delete removes the folder if empty.
            delete@data.FileData(this);
        end
        
        function [x, ctime] = getTrajectory(this, mu, inputidx)
            % Gets a system's trajectory for the given `\mu` and
            % inputindex.
            % Returns [] if no trajectory is found in the Data's Snapshots.
            % 
            % See also: ModelData/getSampleIndex
            
            if nargin < 3
                inputidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            
            x = []; ctime = Inf;
            key = general.Utils.getHash([mu; inputidx]);
            if this.hm.containsKey(key)
                file = this.hm.get(key);
                try
                    s = load(this.getfile(file),'x','ctime');
                    ctime = s.ctime;
                catch%#ok
                    s = load(this.getfile(file),'x');
                    ctime = Inf;
                end
                x = s.x;
            end
        end
        
        function n = getNumTrajectories(this)
           n = this.hm.size;
        end
        
        function [x, mu, inputidx, ctime] = getTrajectoryNr(this, nr)
            % Gets the trajectory with the number nr.
            if nr > this.hm.size || nr < 1
                error('Invalid trajectory number: %d',nr);
            end
            keys = this.hm.keySet.toArray(java_array('java.lang.String',1));
            try
                s = load(this.getfile(this.hm.get(keys(nr))),'x','mu','inputidx','ctime');
                ctime = s.ctime;
                % Workaround for backwards compatibility
            catch%#ok
                s = load(this.getfile(this.hm.get(keys(nr))),'x','mu','inputidx');
                ctime = Inf;
            end
            x = s.x; mu = s.mu; inputidx = s.inputidx;
        end
        
        function addTrajectory(this, x, mu, inputidx, ctime)%#ok
            % Adds a trajectory to the ModelData instance.
            
            if nargin < 4
                inputidx = [];
                if nargin < 3 
                    mu = [];
                end
            end
            
            newd = size(x,1); newmu = length(mu);
            if isempty(this.sizes)
                this.sizes = [newd newmu];
            elseif newd ~= this.sizes(1)
                error('Invalid trajectory dimension. Existing: %d, new: %d',this.sizes(1),newd);
            elseif newmu ~= this.sizes(2)
                error('Invalid parameter dimension. Existing: %d, new: %d',this.sizes(2),newmu);
            end
            
            if isempty(this.trajlen)
                this.trajlen = size(x,2);
            elseif this.trajlen ~= size(x,2)
                error('Invalid trajectory length. Existing: %d, new: %d',this.trajlen,size(x,2));
            end
            
            key = general.Utils.getHash([mu; inputidx]);
            if this.hm.containsKey(key)
                warning('KerMor:MemoryTrajectoryData','Trajectory already present. Replacing.');
            end
            file = [key '.mat'];
            this.hm.put(key,file);
            
            file = fullfile(this.DataDirectory,file);
            try
                save(file,'x','mu','inputidx','ctime');
            catch ME
                this.hm.remove(key);
                rethrow(ME);
            end
            this.updateBB(x);
        end
        
        function clearTrajectories(this)
            ks = this.hm.values.iterator;
            fs = filesep;
            while ks.hasNext
                file = [this.DataDirectory fs ks.next];
                try
                    if exist(file,'file') == 2
                        delete(file);
                    end
                catch ME
                    warning('KerMor:data:FileTrajectoryData','Could not delete file "%s": %s',file,ME.message);
                end
            end
            this.hm.clear;
            this.bbmin = [];
            this.bbmax = [];
            this.sizes = [];
            this.trajlen = [];
        end
        
        function consolidate(this, model, model_ID)
            % Rebuild the hashmap for the current FileData using the current ParamSamples and the models training inputs.
            %
            % This method is used when trajectories are generated within a
            % parfor loop, as then the FileTrajectoryData's are remotely
            % instantiated and their inner hash maps not synced with the
            % local one.
            %
            % Uses the model's parameter samples to compute hashes, look up
            % if the files exist in the given folder corresponding to
            % model_ID (if not given the second arguments ID is used) and
            % re-insert them into the hash map.
            %
            % Parameters:
            % model: The model to use for consolidation (the
            % TrainingInputs property is needed)
            % model_ID: DEBUG PARAM. Allows the consolidation to be
            % performed for a different model ID, for which the model data
            % files have been created. If given, the directory is renamed
            % according to the current model ID.
            if nargin == 3
                olddir = fullfile(KerMor.App.DataStoreDirectory,['rm_' num2str(model_ID)]);
                newdir = fullfile(KerMor.App.DataStoreDirectory,['rm_' num2str(model.ID)]);
                movefile(olddir,newdir);
            end
            this.DataDirectory = fullfile(KerMor.App.DataStoreDirectory,['rm_' num2str(model.ID)]);
            
            this.hm.clear;
            this.bbmin = [];
            this.bbmax = [];
            
            for u=1:max(model.TrainingInputCount,1)
                ui = u;
                if model.TrainingInputCount == 0
                    ui = [];
                end
                for n=1:this.SampleCount
                    mu = this.ParamSamples(:,n);
                    key = general.Utils.getHash([mu; ui]);
                    file = [key '.mat'];
                    ffile = fullfile(this.DataDirectory, file);
                    if exist(ffile,'file') == 2
                        this.hm.put(key,file);
                        % Update the bounding box!
                        s = load(ffile);
                        this.updateBB(s.x);
                    end
                end
            end
        end
        
        function [x,X] = getBoundingBox(this)
            % Gets the bounding box of the state space of all trajectories.
            x = this.bbmin;
            X = this.bbmax;
        end
        
        function [d, mud] = getTrajectoryDoFs(this)
            if isempty(this.sizes)
                if this.getNumTrajectories > 0
                    [x,mu] = this.getTrajectoryNr(1);
                    d = size(x,1);
                    mud = length(mu);
                else
                    d = 0; mud = 0;
                end
                this.sizes = [d mud]; 
            end
            d = this.sizes(1);
            mud = this.sizes(2);
        end
    end
    
    %% general.ABlockSVD implementations
    methods(Access=protected)
        function n = getNumBlocks(this)
            n = this.getNumTrajectories;
        end
        
        function n = getColsPerBlock(this)
            n = this.trajlen;
        end
        
        function B = getBlock(this, nr)
            B = this.getTrajectoryNr(nr);
        end
        
        function [n, m] = getTotalSize(this)
            n = this.getTrajectoryDoFs;
            m = this.getNumBlocks * this.trajlen;
        end
    end
    
    methods(Access=private)
        function updateBB(this, x)
            % Compute current bounding box
            [m,M] = general.Utils.getBoundingBox(x);
            if isempty(this.bbmin)
                this.bbmin = m;
                this.bbmax = M;
            else
                this.bbmin = min(this.bbmin,m);
                this.bbmax = max(this.bbmin,M);
            end
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            % Loads a FileTrajectoryData instance.
            %
            % Ensures that the directory associated with this FileTrajectoryData is existent.
            this = loadobj@data.FileData(this);
            if exist(this.DataDirectory,'dir') ~= 7
                this.hm.clear;
                this.bbmin = [];
                this.bbmax = [];
                warning('FileTrajectoryData:load','This FileTrajectoryData instance was created on host "%s" and the DataDirectory cannot be found on the local machine "%s". Clearing FileTrajectoryData.',this.host,KerMor.getHost);
            end
        end
    end
    
    methods(Static)
        function res = test_FileTrajectoryData
            
            model.ID = 'testModelData';
            m = data.FileTrajectoryData(model);

            T = 10;
            res = true;
            for i=1:T;
                tr(:,:,i) = rand(30,50);%#ok
                p(:,i) = rand(4,1);%#ok
                in(i) = i;%#ok
            end

            % Params only
            for i=1:T;
                m.addTrajectory(tr(:,:,i),p(:,i),[],1);
            end
            res = res && m.getNumTrajectories == T;

            pos = false(1,1,T);
            ppos = false(1,T);
            for i=1:T;
                [x, pi] = m.getTrajectoryNr(i);
                pos = pos | sum(sum(repmat(x,[1 1 T]) - tr,1),2) == 0;
                ppos = ppos | sum(repmat(pi,1,T) - p,1) == 0;
                
                x = m.getTrajectory(p(:,i),[]);
                res = res && isequal(x,tr(:,:,i));
            end
            res = res && all(pos) && all(ppos);
            m.clearTrajectories;

            % Inputs only
            for i=1:T;
                m.addTrajectory(tr(:,:,i),[],in(i),1);
            end
            res = res && m.getNumTrajectories == T;

            pos = false(1,1,T);
            ipos = false(1,T);
            for i=1:T
                [x, ~, ini] = m.getTrajectoryNr(i);
                pos = pos | sum(sum(repmat(x,[1 1 T]) - tr,1),2) == 0;
                ipos = ipos | repmat(ini,1,T)-in == 0;
                
                x = m.getTrajectory([],in(i));
                res = res && isequal(x,tr(:,:,i));
            end
            res = res && all(pos) && all(ipos);
            m.clearTrajectories;
            
            % Both
            for i=1:T;
                m.addTrajectory(tr(:,:,i),p(:,i),in(i),1);
            end
            res = res && m.getNumTrajectories == T;
            
            pos = false(1,1,T);
            ppos = false(1,T);
            ipos = false(1,T);
            for i=1:T
                [x, pi, ini] = m.getTrajectoryNr(i);
                pos = pos | sum(sum(repmat(x,[1 1 T]) - tr,1),2) == 0;
                ppos = ppos | sum(repmat(pi,1,T) - p,1) == 0;
                ipos = ipos | repmat(ini,1,T)-in == 0;
                
                x = m.getTrajectory(p(:,i),in(i));
                res = res && isequal(x,tr(:,:,i));
            end
            res = res && all(pos) && all(ppos) && all(ipos);
            m.clearTrajectories;
        end
    end
    
end