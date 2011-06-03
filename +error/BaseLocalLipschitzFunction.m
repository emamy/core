classdef BaseLocalLipschitzFunction < KerMorObject & ISimConstants & ICloneable
% BaseLocalLipschitzFunction: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @new{0,4,dw,2011-05-31} Implemented/inheriting ISimConstants and ICloneable.
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=protected)
        bellfcn;
    end
    
    methods
        
        function this = BaseLocalLipschitzFunction(bellfunc)
            this = this@KerMorObject;
            this.bellfcn = bellfunc;
        end
        
        % No clone method needed here yet since the bell function is NOT cloned multiple error
        % estimators must react correctly if used for the same model, even though the kernel
        % properties for a reduced model should not change anymore..)
        
    end
    
    methods(Abstract)
        % Evaluates the local lipschitz estimation function.
        %
        % Parameters:
        % di: The distances `d_i(t) = \norm{Vz(t)-x_i}` at time `t`
        % Ct: The coarse error bound `C(t)`. `C\equiv\infty` is allowed, too.
        % t: The current time `t\in[0,T]`
        % mu: The current parameter `\mu`
        ci = evaluate(this,di,Ct,t,mu);
    end
    
end