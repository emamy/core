classdef LocalLipKernelEstimator < error.BaseEstimator
    %COMPWISEKERNELESTIMATOR Summary of this class goes here
    %
    %   Requirements:
    %   ReducedModel.System.f is an instance of AKernelCoreFun
    %   FullModel.Approx is an instance of CompwiseKernelCoreFun
    
    properties(Access=private)
        M1;
        M2;
        M3;
        Ma_norms;
        xi;
    end
    
    properties
        KernelLipschitzFcn;
    end
    
    methods
        function this = LocalLipKernelEstimator(rmodel)
            
            % Validity checks
            msg = error.LocalLipKernelEstimator.validModelForEstimator(rmodel);
            if ~isempty(msg)
                error(msg);
            end
            
            % Call superclass constructor with model argument
            this = this@error.BaseEstimator(rmodel);
            k = rmodel.System.f.SystemKernel;
            
            if isa(k,'kernels.BellFunction')
                this.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;
            else
                this.KernelLipschitzFcn = @k.getGlobalLipschitz;
            end
            this.ExtraODEDims = 3;
            
            fm = rmodel.FullModel;
            % Obtain the correct snapshots
            if ~isempty(fm.Approx)
                % Standard case: the approx function is a kernel expansion.
                
                % Get centers of full approximation
                this.xi = fm.Approx.snData.xi;
                % Get full d x N coeff matrix of approx function
                Ma = fm.Approx.Ma;
            else
                % This is the also possible case that the full core
                % function of the system is a KernelExpansion.
                
                % Get centers of full core function
                this.xi = fm.System.f.snData.xi;
                % Get full d x N coeff matrix of core function
                Ma = fm.System.f.Ma;
            end
            if rmodel.System.f.RotationInvariantKernel
                this.xi = rmodel.System.f.snData.xi;
            end
            
            % Precompute norms
            this.Ma_norms = sqrt(sum(Ma.^2));
            
            % Only prepare matrices if projection was used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                % P = (I-VW^t)
                P = (eye(size(rmodel.V,1)) - rmodel.V*rmodel.W');
                D = P' * rmodel.G * P;
                
                hlp = Ma'*D*Ma;
                this.M1 = hlp;
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(rmodel.FullModel.System.B)
                    try
                        B = rmodel.FullModel.System.B.evaluate([],[]);
                    catch ME%#ok
                        B = rmodel.FullModel.System.B.evaluate(0,rmodel.System.getRandomParam);
                        warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                    end
                    % Filter too small entries
                    hlp = B'*D*B;
                    this.M2 = hlp;
                    hlp = Ma'*D*B;
                    this.M3 = hlp;
                end
            else
                % No projection means no projection error!
                this.M1 = 0;
                this.M2 = 0;
                this.M3 = 0;
            end
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            % extract current error
            eold = x(end-this.ExtraODEDims+1:end);
            e = zeros(3,1);
            
            % Compute \alpha(t)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-this.ExtraODEDims), t, mu);
            % An input function u is set
            if nargin == 5
                e(1) = phi'*this.M1*phi + ut'*this.M2*ut + phi'*this.M3*ut;
                % No input case
            else
                e(1) = phi'*this.M1*phi;
            end
            e(1) = sqrt(max(e(1),0));
            
            % Compute the local lipschitz constant estimations
            if this.ReducedModel.System.f.RotationInvariantKernel
                z = x(1:end-this.ExtraODEDims);
            else
                z = ReducedModel.V*x(1:end-this.ExtraODEDims);
            end
            di = this.xi - repmat(z,1,size(this.xi,2));
            di = sqrt(sum(di.^2,1));
            beta = this.Ma_norms * this.KernelLipschitzFcn(di,Inf,t,mu)';
            e(2) = beta;
            e(3) = eold(1)*beta*exp(-eold(2));
        end
        
        function process(this, t, x, mu, inputidx)%#ok
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('CompWiseErrorEstimator:process','Integral part is all zero. Attention!');
            end
            this.LastError = eint(1,:) + exp(eint(2,:)).*eint(3,:);
        end
        
        function e0 = getE0(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = [this.ReducedModel.getExo(mu); 0; 0];
        end 
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = [];
            if ~isa(rmodel.System.f,'dscomponents.AKernelCoreFun')
                errmsg = 'The reduced model''s core function must be a subclass of dscomponents.AKernelCoreFun for this error estimator.'; 
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

%             ts = tic;
%             beta = this.ReducedModel.FullModel.System.f.getLocalLipschitz(x(1:end-this.ExtraODEDims),t,mu);
%             t = toc(ts);
%             fprintf('beta1: %f, time: %f\n', beta, t);
%             ts = tic;
%             beta2 = this.ReducedModel.FullModel.System.f.getLocalLipschitz2(x(1:end-this.ExtraODEDims),t,mu);
%             t = toc(ts);
%             fprintf('beta2: %f, time: %f\n', beta2, t);
%             ts = tic;
%             beta3 = this.ReducedModel.FullModel.System.f.getLocalLipschitz3(x(1:end-this.ExtraODEDims),t,mu);
%             t = toc(ts);
%             fprintf('beta3: %f, time: %f\n', beta3, t);
%             %beta = min([beta,beta2,beta3]);
%             beta = min([beta,beta2]);
