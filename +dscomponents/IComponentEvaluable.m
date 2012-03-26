classdef IComponentEvaluable < handle
% IComponentEvaluable: 
%
%
%
% @author Daniel Wirtz @date 2012-03-26
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing   
    
    methods(Abstract)
        idx = getComponentArgumentIndices(this, i);
        
        evaluateComponents(this, j, sizes, x, t, mu);
    end
    
    methods
        function res = test_IComponentEvaluableMatch(this, dim, pdim)
            x = rand(dim,1);
            mu = rand(pdim,1);
            t = rand;
            fx = this.evaluateCoreFun(x, t, mu);
            jr = [];
            jend = zeros(dim,1);
            for i=1:dim
                jr = [jr this.getComponentArgumentIndices(i)];%#ok
                jend(i) = length(jr);
            end
            fxc = this.evaluateComponents(1:dim, jend, x(jr), t, mu);
            d = norm(fx-fxc);
            fprintf('Norm difference: %e\n',d);
            res = d < sqrt(eps);
        end
    end
    
end