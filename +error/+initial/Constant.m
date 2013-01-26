classdef Constant < error.initial.Base
% Constant: 
%
%
%
% @author Daniel Wirtz @date 2011-07-04
%
% @new{0,5,dw,2011-07-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        e0;
    end
    
    methods
        function prepareForReducedModel(this, rm)
            fs = rm.FullModel.System;
            x0 = fs.x0.evaluate([]) ./ fs.StateScaling;
            % Only project if projection is used
            if ~isempty(rm.V)
                x0 = x0 - rm.V*(rm.W'*x0);
            end
            this.e0 = Norm.LG(x0,rm.G);
        end
        
        function e0 = getE0(this, mu)%#ok
            e0 = this.e0;
        end
    end
end