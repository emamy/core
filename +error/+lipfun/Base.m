classdef Base < KerMorObject & ICloneable
% Base: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @change{0,5,dw,2011-07-04} Moved all local Lipschitz functions to the package 'error.lipfun'
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
%
% @todo remove old ISimConstants interface implementation
% prepareConstants, at least rename it somehow
    
    properties(Access=protected)
        bellfcn;
    end
    
    methods
        
        function this = Base(bellfunc)
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
        % di: The distances `d_i(t) = \no{Vz(t)-x_i}` at time `t`
        % Ct: The coarse error bound `C(t)`. `C\equiv\infty` is allowed, too.
        ci = evaluate(this, di, Ct);
    end
    
end