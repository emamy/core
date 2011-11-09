classdef MemoryModelData < data.AModelData
% Data class that contains a model's large data, purely in system memory.
%
% @author Daniel Wirtz @date 2011-08-03
%
% @new{0,5,dw,2011-08-03} Added this class (new organization for model data).
%
% @change{0,3,sa,2011-05-10} Implemented setters for the properties
% 
% @change{0,3,dw,2011-04-01} 
% - Changed the old 'ProjTrainData' to 'TrainingData', as this property
% name describes the usage more precisely.
%
% @change{0,1,dw} More common projection via matrices `V,W` instead of
% `V,V^t`.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
%
    
    properties(SetAccess=private)
        % The trajectories in a dim x timelength x trajectorynumber array.
        %
        % Get access allowed for efficient full trajectory data access.
        TrajectoryData;
        
        % The parameters associated with the trajectories. Each column index corresponds to the
        % trajectory number dimension.
        %
        % Empty if none are used.
        Parameters;
        
        % The input indices associated with the trajectories. Each column index corresponds to the
        % trajectory number dimension.
        %
        % Empty if none are used.
        InputIndices;
    end
    
    properties(Access=private)
        % The HashMap used to store the indices for each trajectory.
        hm;
    end
    
    methods
        
        function this = MemoryModelData
            if ~usejava('jvm')
                error('MemoryModelData cannot be used as java is not enabled.');
            end
            this.hm = java.util.HashMap;
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
                idx = this.hm.get(key);
                x = this.TrajectoryData(:,:,idx);
            end
        end
        
        function n = getNumTrajectories(this)
           n = size(this.TrajectoryData,3);
        end
        
        function [x, mu, inputidx] = getTrajectoryNr(this, nr)
            % Gets the trajectory with the number nr.
            if nr > size(this.TrajectoryData,3) || nr < 1
                error('Invalid trajectory number: %d',nr);
            end
            x = this.TrajectoryData(:,:,nr);
            mu = [];
            if ~isempty(this.Parameters)
                mu = this.Parameters(:,nr);
            end
            inputidx = [];
            if ~isempty(this.InputIndices)
                inputidx = this.InputIndices(nr);
            end
        end
        
        function addTrajectory(this, x, mu, inputidx)
            % Adds a trajectory to the ModelData instance.
            
            if nargin < 4
                inputidx = [];
                if nargin < 3 
                    mu = [];
                end
            end
            
            if isempty(inputidx) && ~isempty(this.InputIndices)
                error('MemoryModelData can only hold TrajectoryData of one type. Having InputIndices and not at the same time is not allowed.');
            end
            if isempty(mu) && ~isempty(this.Parameters)
                error('MemoryModelData can only hold TrajectoryData of one type. Having parameters and not at the same time is not allowed.');
            end
            if ~isempty(this.TrajectoryData) && (size(x,1) ~= size(this.TrajectoryData,1) || size(x,2) ~= size(this.TrajectoryData,2))
                error('New trajectory size mismatches the size of already stored TrajectoryData.');
            end
            
            key = general.Utils.getHash([mu; inputidx]);
            if this.hm.containsKey(key)
                warning('KerMor:MemoryModelData','Trajectory already present. Replacing');
                this.TrajectoryData(:,:,this.hm.get(key)) = x;
            else
                this.TrajectoryData(:,:,end+1) = x;
                if ~isempty(mu)
                    this.Parameters(:,end+1) = mu;
                end
                if ~isempty(inputidx)
                    this.InputIndices(end+1) = inputidx;
                end
                this.hm.put(key,size(this.TrajectoryData,3));
            end
        end
        
        function clearTrajectories(this)
            this.TrajectoryData = double.empty(0,0,0);
            this.Parameters = double.empty(0,0);
            this.InputIndices = [];
            this.hm.clear;
        end
        
        function [x,X] = getBoundingBox(this)
            [x, X] = general.Utils.getBoundingBox(this.TrajectoryData(:,:));
        end
    end
    
    methods(Static)
        function res = test_MemoryModelData
            
            m = data.MemoryModelData;
            
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

            for i=1:T;
                [x, pi] = m.getTrajectoryNr(i);
                res = res && isequal(x,tr(:,:,i)) && isequal(p(:,i),pi);
                x = m.getTrajectory(p(:,i),[]);
                res = res && isequal(x,tr(:,:,i));
            end
            m.clearTrajectories;

            % Inputs only
            for i=1:T;
                m.addTrajectory(tr(:,:,i),[],in(i));
            end
            res = res && m.getNumTrajectories == T;

            for i=1:T;
                [x, pi, ini] = m.getTrajectoryNr(i);
                res = res && isequal(x,tr(:,:,i)) && isequal(ini,in(i));
                x = m.getTrajectory([],in(i));
                res = res && isequal(x,tr(:,:,i));
            end
            m.clearTrajectories;

            % both
            for i=1:T;
                m.addTrajectory(tr(:,:,i),p(:,i),in(i));
            end
            res = res && m.getNumTrajectories == T;

            for i=1:T;
                [x, pi, ini] = m.getTrajectoryNr(i);
                res = res && isequal(x,tr(:,:,i)) && isequal(ini,in(i)) && isequal(p(:,i),pi);
                x = m.getTrajectory(p(:,i),in(i));
                res = res && isequal(x,tr(:,:,i));
            end
            m.clearTrajectories;
        end
    end
end

