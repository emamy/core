classdef FileModelData < data.AModelData
% FileModelData: Trajectory datadir in external files.
%
% The constructor takes an optional storage_root parameter.
% If given, it must be a valid folder in the file system.
% If not given, the value set in KerMor.App.DataStoreDirectory is used.
%
% @author Daniel Wirtz @date 2011-08-04
%
% @new{0,5,dw,2011-11-02} Implemented the getBoundingBox method from superclass. Now keeping
% track of bounding box while adding trajectories.
%
% @new{0,5,dw,2011-10-14} Added a new consolidate method in order to
% rebuild the internal index hashmap from the files in a directory and the
% data.AModelData.ParamSamples. This method is called after parallel
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
    
    properties%(Access=private)
        % The HashMap used to store the indices for each trajectory.
        hm;
        
        datadir;
    end
    
    properties(Access=private)
        % Bounding box minimum
        bbmin = [];
        
        % Bounding box maximum
        bbmax = [];
    end
    
    methods
        function this = FileModelData(model, storage_root)
            % Creates a new ModelData instance with trajectory datadir in a file folder.
            %
            % Parameters:
            % model: The model the data is stored for.
            % storage_root: [Optional] If given, it must be a valid folder in the file system. If
            % not given, the value set in KerMor.App.DataStoreDirectory is used. @type char
            if ~usejava('jvm')
                error('FileModelData cannot be used as java is not enabled.');
            end
            this.hm = java.util.HashMap;
            if nargin == 2
                if isa(storage_root,'char') && exist(storage_root,'dir') == 7
                    this.datadir = storage_root;
                else
                    error('Invalid folder: %s',storage_root);
                end
            else
                this.datadir = KerMor.App.DataStoreDirectory;
            end
            this.datadir = fullfile(this.datadir,['rm_' num2str(model.ID)]);
            if exist(this.datadir,'dir') ~= 7
                try
                    mkdir(this.datadir);
                catch ME
                    me = MException('KerMor:data:FileModelData','Could not create dir "%s"',this.datadir);
                    me.addCause(ME);
                    me.throw;
                end
            end
            this.clearTrajectories;
        end
        
        function x = getTrajectory(this, mu, inputidx)
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
            
            x = [];
            key = general.Utils.getHash([mu; inputidx]);
            if this.hm.containsKey(key)
                file = this.hm.get(key);
                s = load(this.getfile(file),'x');
                x = s.x;
            end
        end
        
        function n = getNumTrajectories(this)
           n = this.hm.size;
        end
        
        function [x, mu, inputidx] = getTrajectoryNr(this, nr)
            % Gets the trajectory with the number nr.
            if nr > this.hm.size || nr < 1
                error('Invalid trajectory number: %d',nr);
            end
            keys = this.hm.keySet.toArray(java_array('java.lang.String',1));
            s = load(this.getfile(this.hm.get(keys(nr))),'x','mu','inputidx');
            x = s.x; mu = s.mu; inputidx = s.inputidx;
        end
        
        function addTrajectory(this, x, mu, inputidx)
            % Adds a trajectory to the ModelData instance.
            
            if nargin < 4
                inputidx = [];
                if nargin < 3 
                    mu = [];
                end
            end
            
            key = general.Utils.getHash([mu; inputidx]);
            if this.hm.containsKey(key)
                warning('KerMor:MemoryModelData','Trajectory already present. Replacing.');
            end
            file = [key '.mat'];
            this.hm.put(key,file);
            
            file = fullfile(this.datadir,file);
            try
                save(file,'x','mu','inputidx');
            catch ME
                this.hm.remove(key);
                rethrow(ME);
            end
            this.updateBB(x);
        end
        
        function clearTrajectories(this)
            ks = this.hm.values.iterator;
            while ks.hasNext
                file = fullfile(this.datadir,ks.next);
                try
                    delete(file);
                catch ME
                    warning('KerMor:data:FileModelData','Could not delete file "%s": %s',file,ME.message);
                end
            end
            this.hm.clear;
            this.bbmin = [];
            this.bbmax = [];
        end
        
        function consolidate(this, model, model_ID)
            % Rebuild the hashmap for the current FileData using the current ParamSamples and the models training inputs.
            %
            % This method is used when trajectories are generated within a
            % parfor loop, as then the FileModelData's are remotely
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
            this.datadir = fullfile(KerMor.App.DataStoreDirectory,['rm_' num2str(model.ID)]);
            
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
                    ffile = fullfile(this.datadir, file);
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
    end
    
    methods(Access=private)
        function file = getfile(this, file)
            file = fullfile(this.datadir,file);
            if exist(file,'file') ~= 2
                error('File not found: "%s". Have you deleted model data files?',file);
            end
        end
        
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
    
    methods(Static)
        function res = test_FileModelData
            
            model.ID = 'testModelData';
            m = data.FileModelData(model);

            T = 10;
            res = true;
            for i=1:T;
                tr(:,:,i) = rand(30,50);%#ok
                p(:,i) = rand(4,1);%#ok
                in(i) = i;%#ok
            end

            % Params only
            for i=1:T;
                m.addTrajectory(tr(:,:,i),p(:,i),[]);
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
                m.addTrajectory(tr(:,:,i),[],in(i));
            end
            res = res && m.getNumTrajectories == T;

            pos = false(1,1,T);
            ipos = false(1,T);
            for i=1:T
                [x, pi, ini] = m.getTrajectoryNr(i);
                pos = pos | sum(sum(repmat(x,[1 1 T]) - tr,1),2) == 0;
                ipos = ipos | repmat(ini,1,T)-in == 0;
                
                x = m.getTrajectory([],in(i));
                res = res && isequal(x,tr(:,:,i));
            end
            res = res && all(pos) && all(ipos);
            m.clearTrajectories;
            
            % Both
            for i=1:T;
                m.addTrajectory(tr(:,:,i),p(:,i),in(i));
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