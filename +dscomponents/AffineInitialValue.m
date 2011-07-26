classdef AffineInitialValue < dscomponents.AInitialValue & general.AffParamMatrix
% AffineInitialValue: Parameter-affine initial value for dynamical systems.
%
% Extends the standard AffParamMatrix from the general package.
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
    
    methods
        function x0 = evaluate(this, mu)
            x0 = this.compose(mu);
        end
        
        function x0 = compose(this, mu)
            % Call affine function evaluate with empty time
            x0 = compose@general.AffParamMatrix(this, [], mu);
        end
        
        function projected = project(this, V, W)%#ok
            % Projects the affine parametric initial value into the subspace spanned by
            % `V,W`.
            
            % Uses the overridden operators in AffParamMatrix to create a copy.
            projected = W'*this;
        end
        
        function copy = clone(this)
            copy = dscomponents.AffineInitialValue;
            copy = clone@general.AffParamMatrix(this, copy);
        end
    end
    
end