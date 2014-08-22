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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        function [xi, ti, mui, fxi] = select(this, model)
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
            % fxi: The `f(x_i,t_i,\mu_i)` evaluations. Empty if not
            % computed yet. @type matrix
            
            md = model.Data;
            td = md.TrajectoryData;
            hasfxi = false;
            if ~isempty(md.TrajectoryFxiData) && td.getNumTrajectories == md.TrajectoryFxiData.getNumTrajectories
                hasfxi = true;
            end
            nt = td.getNumTrajectories;
            
            % The trajectory length
            tl = length(model.Times);
            % The total size of all trajectories
            % If nonuniform data is allowed, we need to iterate over all
            % data in order to determine the different lengths
            if isa(td,'data.FileTrajectoryData') && ~td.UniformTrajectories
                sizes = zeros(1,nt);
                pi = ProcessIndicator('Selecting %d approximation training data samples from %d trajectories',nt,false,this.Size,nt);
                for n = 1:nt
                    x = td.getTrajectoryNr(n);
                    sizes(n) = size(x,2);
                end
            else
                sizes = ones(1,nt) * tl;
            end
            
            linpos = [1 cumsum(sizes)];
            ts = linpos(end);
            
            % Select maximal as many as are there!
            s = min(this.Size,ts);
            idx = round(linspace(1,ts,s));
            
            [xd, mud] = td.getTrajectoryDoFs;
            % Use 512 MB chunks for approx train data
            xi = data.FileMatrix(xd,s,'Dir',md.DataDirectory);
            ti = zeros(1,s);
            mui = zeros(mud,s);
            fxi = [];
            if hasfxi
                fxi = data.FileMatrix(xd,s,'Dir',md.DataDirectory);
            end
            
            pi = ProcessIndicator('Selecting %d approximation training data samples from %d trajectories',nt,false,this.Size,nt);
            atdpos = 0;
            for k=1:nt
                sel = idx(idx >= linpos(k) & idx < linpos(k+1));
                if ~isempty(sel)
                    sel = sel - min(sel) + 1;
                    [x, mu] = td.getTrajectoryNr(k);
                    atdpos = atdpos(end) + (1:length(sel));
                    xi(:,atdpos) = x(:,sel);
                    ti(atdpos) = model.Times(sel);
                    mui(:,atdpos) = repmat(mu,1,length(atdpos));

                    if hasfxi
                        [fx, fxmu] = md.TrajectoryFxiData.getTrajectoryNr(k);
                        if ~isequal(mu,fxmu)
                            hasfxi = false;
                        else
                            fxi(:,atdpos) = fx(:,sel);%#ok
                        end
                    end
                end
                pi.step;
            end
            pi.stop;
        end
    end
    
    methods(Static)
        function res = test_LinSpaceSelector
            [res, m] = models.burgers.Tests.test_Burgers_DEIM_versions(50,2);
            m.Approx.TrainDataSelector = data.selection.LinspaceSelector;
            m.Approx.TrainDataSelector.Size = 300;
            m.off4_genApproximationTrainData;
            m.ComputeTrajectoryFxiData = true;
            m.off2_genTrainingData;
            m.off4_genApproximationTrainData;
            res = res & true;
        end
    end
end

