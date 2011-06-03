classdef GLEstimator < error.BaseKernelEstimator
    % GLEstimator: Global lipschitz constant error estimator
    %
    % @author Daniel Wirtz @date 2010-05-10
    %
    % @change{0,4,dw,2011-05-29} Restructured the error estimators to better adopt to the current
    % formulation. Now the KernelEstimators have a function getBeta instead of implementing the
    % evalODEPart by themselves.
    %
    % @change{0,4,dw,2011-05-25} Changed the implementation to correspond to the new comparison
    % lemma estimator derivation. This requires only one ODE dimension and is a sharper estimate
    % than before.
    %
    % @change{0,4,dw,2011-05-23} Adopted to the new error.BaseEstimator interface with separate output
    % error computation.
    
    properties(Access=private)
        cf;
    end
    
    methods
        function this = GLEstimator(rmodel)
            this = this@error.BaseKernelEstimator;
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.GLEstimator;
            copy = clone@error.BaseKernelEstimator(this, copy);
            copy.cf = this.cf;
        end
        
        function prepareConstants(this)
            prepareConstants@error.BaseKernelEstimator(this);
            % Standard case: the approx function is a kernel expansion. it
            % can also be that the system's core function is already a
            % kernel expansion
            fm = this.ReducedModel.FullModel;
            if ~isempty(fm.Approx)
                % Get full d x N coeff matrix of approx function
                f = fm.Approx;
            else
                % Get full d x N coeff matrix of core function
                f = fm.System.f;
            end
            this.cf = f.getGlobalLipschitz(0, []);
        end
    end
    
    methods(Access=protected)
        function b = getBeta(this, x, t, mu)%#ok
            b = this.cf;
        end
        
        function postprocess(this, t, x, mu, inputidx)
            postprocess@error.BaseKernelEstimator(this, t, x, mu, inputidx);
            this.StateError(1,:) = x(end,:);
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = validModelForEstimator@error.BaseKernelEstimator(rmodel);
            if isempty(errmsg) && ~isa(rmodel.FullModel.System.f,'dscomponents.IGlobalLipschitz')
                errmsg = 'The full model''s core function must implement the dscomponents.IGlobalLipschitz interface for this error estimator.';
            end
        end
    end
    
end