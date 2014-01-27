classdef IReductionSummaryPlotProvider < handle
% IReductionSummaryPlotProvider: 
%
%
%
% @author Daniel Wirtz @date 2012-09-25
%
% @new{0,6,dw,2012-09-25} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods(Abstract)
        plotSummary(this, pm, context);
    end
    
end