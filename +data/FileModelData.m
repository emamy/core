classdef FileModelData < data.AModelData
% FileModelData: Trajectory datadir in external files.
%
% The constructor takes an optional storage_root parameter.
% If given, it must be a valid folder in the file system.
% If not given, the value set in KerMor.App.DataStoreDirectory is used.
%
% @author Daniel Wirtz @date 2011-08-04
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
        
        datadir;
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
        
        function addTrajectory(this, x, mu, inputidx)%#ok
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
        end
    end
    
    methods(Access=private)
        function file = getfile(this, file)
            file = fullfile(this.datadir,file);
            if exist(file,'file') ~= 2
                error('File not found: "%s". Have you deleted model data files?',file);
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