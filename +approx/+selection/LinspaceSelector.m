classdef LinspaceSelector < approx.selection.ASelector
% Selects Size equally spaced samples from the training data.
% 
% Equally spaced in this context means not spaced in a state space sense.
% Instead, the index set resulting from concatenation of ALL trajectories available will be sampled
% over Size linearly equidistant indices. This way, all trajectories are treated as if it was one
% big one.
%
% @author Daniel Wirtz @date 2011-04-12
%
% @change{0,5,dw,2011-08-4} Changed the behaviour of this class to adopt to the new AModelData
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
    end
    
    methods
        
        function this = LinspaceSelector
            this.registerProps('Size');
        end
        
        function copy = clone(this)
            copy = approx.selection.LinspaceSelector;
            copy = clone@approx.selection.ASelector(this, copy);
            copy.Size = this.Size;
        end
        
        function set.Size(this, value)
            if ~isposintscalar(value)
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
            
            nt = model.Data.getNumTrajectories;
            
            % The trajectory length
            tl = length(model.Times);
            % The total size of all trajectories
            ts = nt*tl;
            % Get valid size
            s = min(this.Size,ts);
            % Obtain indices for each trajectory
            idx = mod(round(linspace(0,ts-1,s)),tl)+1;
            pos = [0 find(idx(2:end)-idx(1:end-1) < 0) s];
            
            xi = [];
            ti = [];
            mui = [];
            for k=1:nt
                [x, mu] = model.Data.getTrajectoryNr(k);
                sel = pos(k)+1:pos(k+1);
                xi = [xi x(:,idx(sel))]; %#ok<*AGROW>
                ti = [ti model.Times(idx(sel))];
                mui = [mui repmat(mu,1,length(sel))];
            end
        end
    end
end

