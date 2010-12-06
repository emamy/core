classdef GlobalLipKernelEstimator < error.BaseEstimator
    %COMPWISEKERNELESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function this = GlobalLipKernelEstimator(rmodel)
            % @todo: check if the computations for M1,M2 etc can already be
            % done at the time of projection?
            
            this.ExtraODEDims = 1;
            
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            
            copy = error.GlobalLipKernelEstimator;
            % copy.KernelLipschitzFcn/ODEDims is assigned in constructor!
            copy = clone@error.BaseEstimator(this, copy);
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-this.ExtraODEDims,:), t, mu);
            % An input function u is set
            if nargin == 5
                e = phi'*this.M1*phi + ut'*this.M2*ut + phi'*this.M3*ut;
                % No input case
            else
                e = phi'*this.M1*phi;
            end
            e = sqrt(max(e,0));
        end
        
        function process(this, t, x, mu, inputidx)%#ok
            % @todo also take into account FullModel.Approx at cf computation!
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('CompWiseErrorEstimator:process','Integral part is all zero. Attention!');
            end
            cf = this.ReducedModel.FullModel.System.f.getGlobalLipschitz(t,mu);
            pt1 = exp(cf .* reshape(t,1,[]));
            pt2 = (eint + this.ReducedModel.getExo(mu));
            this.LastError = pt1 .* pt2;
        end
        
        function e0 = getE0(this, mu)%#ok
            % Returns the initial error at `t=0` of the integral part.
            % For this estimator, this is simply one dimension and zero.
            e0 = 0;
        end
        
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = [];
            if ~isa(rmodel.FullModel.System.f,'dscomponents.IGlobalLipschitz')
                errmsg = 'The full model''s core function must implement the dscomponents.IGlobalLipschitz interface for this error estimator.';
            elseif ~isempty(rmodel.FullModel.Approx) 
                if ~isa(rmodel.FullModel.Approx,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'The full model''s approx function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
                end
            elseif ~isa(rmodel.FullModel.System.f,'dscomponents.CompwiseKernelCoreFun')
                errmsg = 'If no approximation is used, the full model''s core function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
            end
        end
    end
    
end