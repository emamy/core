classdef GLEstimator < error.BaseCompLemmaEstimator
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
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing    
    %
    % @todo extend the GlobalLipschitz computation to also depend on t,mu (with check ahead if
    % constant, then efficient one-time eval is possible)
    
    properties(Access=private)
        cf;
        f;
    end
    
    methods
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.GLEstimator;
            copy = clone@error.BaseCompLemmaEstimator(this, copy);
            copy.cf = this.cf;
            copy.f = this.f;
        end
        
        function b = getBeta(this, x, t)%#ok
            if ~isempty(this.cf)
                b = this.cf;
            else
                b = this.f.getGlobalLipschitz(t, this.mu);
            end
        end
        
        function ct = prepareConstants(this, mu, inputidx)
            % Return values:
            % ct: The time needed for preprocessing @type double
            
            ct = prepareConstants@error.BaseCompLemmaEstimator(this, mu, inputidx);
            % Standard case: the approx function is a kernel expansion. it
            % can also be that the system's core function is already a
            % kernel expansion
            st = tic;
            fm = this.ReducedModel.FullModel;
            if ~isempty(fm.Approx)
                % Get full d x N coeff matrix of approx function
                this.f = fm.Approx;
            else
                % Get full d x N coeff matrix of core function
                this.f = fm.System.f;
            end
            % Only precompute the global lipschitz constant if it is (t,\mu)-independent
            if isa(this.f,'kernels.KernelExpansion') || (isa(this.f,'kernels.ParamTimeKernelExpansion') && ...
                    isa(this.f.ParamKernel,'kernels.NoKernel') && isa(this.f.TimeKernel,'kernels.NoKernel'))
                this.cf = this.f.getGlobalLipschitz(0, mu);
            else
                this.cf = [];
            end
            ct = ct + toc(st);
        end
        
        function ct = postProcess(this, x, t, inputidx)
            % Return values:
            % ct: The time needed for postprocessing @type double
            
            st = tic;
            this.StateError(1,:) = x(end,:);
            ct = postProcess@error.BaseCompLemmaEstimator(this, t, x, mu, inputidx) + toc(st);
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(model)
            % Validations
            errmsg = validModelForEstimator@error.BaseCompLemmaEstimator(model);
            if isempty(errmsg) && ~isa(model.System.f.Expansion,'dscomponents.IGlobalLipschitz')
                errmsg = 'The full model''s core function must implement the dscomponents.IGlobalLipschitz interface for this error estimator.';
            end
        end
    end
    
end