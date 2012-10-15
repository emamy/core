classdef BaseSpaceReducer < KerMorObject & general.IReductionSummaryPlotProvider
    % Base class for all space reduction algorithms.
    %
    % @author Daniel Wirtz
    % @date 11.03.2010
    
    properties%(SetAccess=protected)
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
    
end

