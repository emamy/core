classdef BaseSpaceReducer < KerMorObject & IReductionSummaryPlotProvider
    % Base class for all space reduction algorithms.
    %
    % @author Daniel Wirtz
    % @date 11.03.2010
    
    properties(SetObservable)
        % Flag if to include `f(x_i)` values for each `x_i` value, too
        %
        % @propclass{important} Including this data in the subspace can reduce projection
        % errors of nonlinearities
        %
        % @type logical @default false
        IncludeTrajectoryFxiData = false;
        
        IncludeFiniteDifferences = false;
        
        IncludeBSpan = false;
        
        IncludeAxData = false;
        
        % Set this to an index vector to have the space reducer
        % reduce only the indicated dimensions.
        %
        % The resulting reduced space will at least be as large as the
        % pertained dimensions and V (W) will contain unit vectors at the
        % excluded dimensions. Set to ':' to get full reduction.
        %
        % @type colvec<integer> @default ':'
        ReducableDims = ':';
    end
    
    properties(SetAccess=protected)
        ProjectionError = [];
    end
    
    methods
        
        function [V, W] = generateReducedSpace(this, model)
            [V,W] = this.generateReducedSpaceImpl(model);
            % For partial reduction we need to insert identity for the
            % pertained dimensions.
            if ~strcmp(this.ReducableDims,':')
                dim = model.System.StateSpaceDimension;
                V = augment(V);
                if ~isempty(W)
                    W = augment(W);
                end
            end
            
            function Afull = augment(A)
                basesize = size(A,2);
                Afull = zeros(dim,basesize);
                Afull(this.ReducableDims,:) = A;
                
                unreduced = 1:dim;
                unreduced(this.ReducableDims) = [];
                nunred = length(unreduced);
                Afull(unreduced,basesize+1:basesize+nunred) = eye(nunred);
            end
        end
        
        function plotSummary(this, pm, context)
            if ~isempty(this.ProjectionError)
                str = sprintf('%s/%s: Projection error over training data',...
                    context,class(this));
                h = pm.nextPlot('spacereduction_projerr',str,'subspace size','error');
                semilogy(h,this.ProjectionError,'LineWidth',2);
            end
        end
    end
    
    methods(Access=protected)
        
        function V = getInitialSpace(this, blockdata, pod, reducable)
            % Computes the initial space, which is the first POD mode of
            % the initial values!
            
            n = blockdata.getNumBlocks;
            x = blockdata.getBlock(1);
            x0 = x(reducable,1);
            if this.ComputeParallel
                parfor idx=2:n
                    x = blockdata.getBlock(idx);%#ok
                    x0 = [x0, x(reducable,1)];
                end
                x0 = unique(x0','rows')';
            else
                for idx=2:n
                    x = blockdata.getBlock(idx);
                    x = x(reducable,1);
                    % Only add nonexisting vectors
                    if isempty(Utils.findVecInMatrix(x0,x))
                        x0 = [x0 x];%#ok
                    end
                end
            end
            if all(x0 == 0)
                if KerMor.App.Verbose > 1
                    fprintf('Initial values are all zero vectors. Using main POD mode of first block data as initial space.\n');
                end
                xb1 = blockdata.getBlock(1);
                V = pod.computePOD(xb1(reducable,:));
            elseif size(x0,2) > 1
                V = pod.computePOD(x0);
            else
                V = x0;
            end
            V = V / norm(V);
        end
    end
        
    methods(Access=protected, Abstract)
        [V,W] = generateReducedSpaceImpl(this, model);
    end
    
    %% Getter & setter
    methods 
        function set.IncludeTrajectoryFxiData(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeTrajectoryFxiData = value;
        end
        
        function set.IncludeFiniteDifferences(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeFiniteDifferences = value;
        end
        
        function set.IncludeBSpan(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeBSpan = value;
        end
        
        function set.IncludeAxData(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeAxData = value;
        end
    end
    
end

