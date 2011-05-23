classdef LocalSecantLipschitz < error.BaseLocalLipschitzFunction
% LocalSecantLipschitz: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = LocalSecantLipschitz(bellfcn)
            this = this@error.BaseLocalLipschitzFunction(bellfcn);
        end
        
        function ci = evaluate(this, di, C, t, mu)%#ok
            b = this.bellfcn;
            x0 = b.x0;
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                ci = ones(size(di))*abs(b.evaluateD1(x0));
            else
                case1 = di - C - x0 > 0;
                case2 = di + C - x0 < 0;
                
                % If C is too small we just take the derivative at this
                % point.
                if C < sqrt(eps)
                    ci(case1 | case2) = abs(b.evaluateD1(di(case1 | case2)));
                else
                    ci(case1) = (b.evaluateScalar(di(case1)-C) - b.evaluateScalar(di(case1))) / C;
                    ci(case2) = (b.evaluateScalar(di(case2)) - b.evaluateScalar(di(case2)+C)) / C;
                end
                ci(~case1 & ~case2) = abs(b.evaluateD1(x0));
            end
        end
    end
    
end