classdef BaseSpaceReducer < KerMorObject & general.IReductionSummaryPlotProvider
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
    end
    
    properties(SetAccess=protected)
        ProjectionError = [];
    end
    
    methods
        function plotSummary(this, pm, context)
            if ~isempty(this.ProjectionError)
                str = sprintf('%s/%s: Projection error over training data',...
                    context,class(this));
                h = pm.nextPlot('spacereduction_projerr',str,'subspace size','error');
                semilogy(h,this.ProjectionError,'LineWidth',2);
            end
        end
    end
        
    methods(Abstract)
        [V,W] = generateReducedSpace(model);
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

