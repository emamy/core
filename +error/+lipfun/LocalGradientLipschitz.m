classdef LocalGradientLipschitz < error.lipfun.Base
% LocalGradientLipschitz: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @new{0,4,dw,2011-05-31} Added new prepareConstants init function from Base.
%
% @change{0,4,dw,2011-05-29} Changed `x_0` to `r_0` to adopt new notation.
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    properties(Access=private,Transient)
        d1;
    end
    
    methods
        function this = LocalGradientLipschitz(bellfcn)
            this = this@error.lipfun.Base(bellfcn);
        end
        
        function copy = clone(this)
            copy = error.lipfun.LocalSecantLipschitz(this.bellfcn);
            %copy = clone@error.lipfun.Base(this, copy);
            copy.d1 = this.d1;
        end
        
        function prepareConstants(this)
            this.d1 = abs(this.bellfcn.evaluateD1(r0));
        end
        
        function ci = evaluate(this, di, C)
            % Consider C=Inf case separately for speed reasons
            b = this.bellfcn;
            r0 = b.r0;
            if isinf(C)
                ci = ones(size(di))*this.d1;
            else
                case1 = di - C - r0 > 0;
                case2 = di + C - r0 < 0;
                ci(case1) = abs(b.evaluateD1(di(case1)-C));
                ci(case2) = abs(b.evaluateD1(di(case2)+C));
                ci(~case1 & ~case2) = this.d1;
            end
        end
    end 
end