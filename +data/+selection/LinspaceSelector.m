classdef LinspaceSelector < data.selection.ASelector
% Selects Size equally spaced samples from the training data.
% 
% Equally spaced in this context means not spaced in a state space sense.
% Instead, the index set resulting from concatenation of ALL trajectories available will be sampled
% over Size linearly equidistant indices. This way, all trajectories are treated as if it was one
% big one.
%
% @author Daniel Wirtz @date 2011-04-12
%
% @change{0,5,dw,2011-08-4} Changed the behaviour of this class to adopt to the new ATrajectoryData
% iterface. The sampling is done as before, treating all trajectories as if they were contained in
% one big array.
%
% @new{0,4,dw,2011-05-06}
% - Added this class to @ref propclasses.
% - Implemented ICloneable interface.
%
% @new{0,3,dw,2011-04-12} Added this class to allow for homogeneous training data selection.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The (maximum) number of elements to take
        %
        % @propclass{critical} The amount of approximation training data to take.
        %
        % @default `\infty`
        Size = Inf;
        
        MaxMemPerc = .4;
    end
    
    methods
        
        function this = LinspaceSelector
            this.registerProps('Size');
        end
        
        function copy = clone(this)
            copy = data.selection.LinspaceSelector;
            %copy = clone@data.selection.ASelector(this, copy);
            copy.Size = this.Size;
        end
        
        function set.Size(this, value)
            if ~isposintscalar(value) && ~isinf(value)
                error('The value must be a positive integer.');
            end
            this.Size = value;
        end
    end
    
    methods(Access=protected,Sealed)
        function [xi, ti, mui] = select(this, model)
            % Selects Size equally (index-)spaced samples from the training data.
            % 
            % Equally spaced in this context means with respect to the indices over the WHOLE
            % trajectories summarized length, NOT to the spatial distances between the training
            % samples.
            %
            % Parameters:
            % model: The full model with the training data @type models.BaseFullModel
            %
            % Return values:
            % xi: The selected `x_i = x(t_i)` training data @type matrix
            % ti: The selected training times `t_i` @type rowvec
            % mui: The selected parameter samples `\mu_i` with which the states
            % `x_i` have been reached @type matrix
            
            md = model.Data;
            td = md.TrajectoryData;
            nt = td.getNumTrajectories;
            
            % The trajectory length
            tl = length(model.Times);
            % The total size of all trajectories
            ts = nt*tl;
            % Get valid size
            s = min(this.Size,ts);
            
            % Obtain indices for each trajectory
            idx = mod(round(linspace(0,ts-1,s)),tl)+1;
            tidx = ceil(linspace(1,ts-1,s)/tl);
            traj = unique(tidx);
            
            [xd, mud] = td.getTrajectoryDoFs;
            xi = data.FileMatrix(xd,s,md.DataDirectory,512*1024^2); % Use 512 MB chunks for approx train data
            ti = zeros(1,s);
            mui = zeros(mud,s);
            
            pi = tools.ProcessIndicator('Selecting approximation training data from %d trajectories',nt,false,nt);
            atdpos = 0;
            for k=1:length(traj)
                [x, mu] = md.TrajectoryData.getTrajectoryNr(traj(k));
                sel = idx(tidx == traj(k));
                atdpos = atdpos(end) + (1:length(sel));
                xi(:,atdpos) = x(:,sel);
                ti(atdpos) = model.Times(sel);
                mui(:,atdpos) = repmat(mu,1,length(atdpos));
                pi.step;
            end
            pi.stop;
        end
    end
end

