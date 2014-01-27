classdef LocalSecantLipschitz < error.lipfun.Base
% LocalSecantLipschitz: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @new{0,4,dw,2011-05-31} Added new init function from Base.
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    properties(Transient, Access=private)
        dx0;
    end
    
    methods
        function this = LocalSecantLipschitz(bellfcn)
            this = this@error.lipfun.Base(bellfcn);
        end
        
        function copy = clone(this)
            copy = error.lipfun.LocalSecantLipschitz(this.bellfcn);
            %copy = clone@error.lipfun.Base(this, copy);
            copy.dx0 = this.dx0;
        end
        
        function prepareConstants(this)
            this.dx0 = abs(this.bellfcn.evaluateD1(x0));
        end
        
        function ci = evaluate(this, di, C)
            b = this.bellfcn;
            x0 = b.x0;
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                ci = ones(size(di))*this.dx0;
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
                ci(~case1 & ~case2) = this.dx0;
            end
        end
    end
    
end