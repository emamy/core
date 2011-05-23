classdef GlobalLipKernelEstimator < error.BaseLipKernelEstimator
    % GlobalLipKernelEstimator: Global lipschitz constant error estimator
    %
    % @author Daniel Wirtz @date 2010-05-10
    %
    % @change{0,4,dw,2011-05-23} Adopted to the new error.BaseEstimator interface with separate output
    % error computation.
    
    methods
        function this = GlobalLipKernelEstimator(rmodel)
            this.ExtraODEDims = 1;
            
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.GlobalLipKernelEstimator;
            copy = clone@error.BaseLipKernelEstimator(this, copy);
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-this.ExtraODEDims,:), t, mu);
            
            % An input function u is set
            if nargin == 5
                e = phi*this.M1*phi' + phi*this.M2*ut + ut'*this.M3*ut;
                % No input case
            else
                e = phi*this.M1*phi';
            end
            e = sqrt(max(e,0));
        end
                
        function e0 = getE0(this, mu)%#ok
            % Returns the initial error at `t=0` of the integral part.
            % For this estimator, this is simply one dimension and zero.
            e0 = 0;
        end
    end
    
    methods(Access=protected)
        function postprocess(this, t, x, mu, inputidx)%#ok
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('CompWiseErrorEstimator:process','Integral part is all zero. Attention!');
            end
            cf = this.ReducedModel.FullModel.System.f.getGlobalLipschitz(t,mu);
            
            this.StateError = exp(cf .* reshape(t,1,[])) .* (eint + this.ReducedModel.getExo(mu));
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = validModelForEstimator@error.BaseLipKernelEstimator(rmodel);
            if isempty(errmsg) && ~isa(rmodel.FullModel.System.f,'dscomponents.IGlobalLipschitz')
                errmsg = 'The full model''s core function must implement the dscomponents.IGlobalLipschitz interface for this error estimator.';
            end
        end
    end
    
end