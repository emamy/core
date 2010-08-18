classdef GlobalLipKernelEstimator < error.BaseEstimator
    %COMPWISEKERNELESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=private)
        M1;
        M2;
        M3;
    end
    
    methods
        function this = GlobalLipKernelEstimator(rmodel)
            % @todo: check if the computations for M1,M2 etc can already be
            % done at the time of projection?
            
            % Validity checks
            msg = error.GlobalLipKernelEstimator.validModelForEstimator(rmodel);
            if ~isempty(msg)
                error(msg);
            end
            
            % Call superclass constructor with model argument
            this = this@error.BaseEstimator(rmodel);
            this.ExtraODEDims = 1;
            
            % Only prepare matrices if projection was used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                % P = (I-VW^t)
                P = (eye(size(rmodel.V,1)) - rmodel.V*rmodel.W');
                D = P' * rmodel.G * P;
                
                fm = rmodel.FullModel;
                % Obtain the correct snapshots
                if ~isempty(fm.Approx)
                    % Standard case: the approx function is a kernel expansion.

                    % Get full d x N coeff matrix of approx function
                    Ma = fm.Approx.Ma;
                else
                    % This is the also possible case that the full core
                    % function of the system is a KernelExpansion.

                    % Get full d x N coeff matrix of core function
                    Ma = fm.System.f.Ma;
                end
                
                this.M1 = Ma'*D*Ma;
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(rmodel.FullModel.System.B)
                    try
                        B = rmodel.FullModel.System.B.evaluate([],[]);
                    catch ME
                        B = rmodel.FullModel.System.B.evaluate(0,rmodel.System.getRandomParam);
                        warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                    end
                    this.M2 = B'*D*B;
                    this.M3 = Ma'*D*B;
                end
                
            else
                % No projection means no error!
                this.M1 = 0;
                this.M2 = 0;
                this.M3 = 0;
            end
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