classdef ExpensiveBetaEstimator < error.BaseLipKernelEstimator
% ExpensiveBetaEstimator: 
%
% This error estimator performs the same estimation process than the local lipschitz estimators,
% however, it computes the full trajectory and makes use of it in the beta computation.
% 
% Two strategies are implemented, one simply computing `\beta(t) = \frac{\no{\fa(x(t)) -
% \fa(Vz(t))}}{\no{x(t)-Vz(t)}}`, so that `\beta(t)` is always the best possible estimation.
% The other one computes the maximum of the two jacobians at `x(t)` and `Vz(t)`. This is not
% rigorous but an estimate (rigorous would be the maximum norm over the whole line from `x(t)` to
% `Vz(t)`)
%
% @author Daniel Wirtz @date 2011-05-18
%
% @new{0,4,dw,2011-05-18} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % Which version to use
        % 1: Direct `beta(t)` computation as explained above
        % 2: Jacobians at `x(t)` and `Vz(t)`
        Version = 1;
    end
    
    properties(Access=private,Transient)
        neg_e1=false;
        xfull;
        fax;
        cnt;
    end
    
    methods
        function this = ExpensiveBetaEstimator(rmodel)
            this.ExtraODEDims = 2;
            if nargin == 1
                this.setReducedModel(rmodel);
            end
            this.cnt = 1;
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.ExpensiveBetaEstimator;
            % ExtraODEDims is set in constructor!
            copy = clone@error.BaseLipKernelEstimator(this, copy);            
            copy.neg_e1 = this.neg_e1;
        end
        
        function clear(this)
            clear@error.BaseLipKernelEstimator(this);
            this.cnt = 1;
            this.neg_e1 = false;
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            % extract current error
            eold = x(end-this.ExtraODEDims+1:end);
            e = zeros(this.ExtraODEDims,1);
            
            % Compute \alpha(t)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-this.ExtraODEDims), t, mu);
            
            e(1) = phi*this.M1*phi';
            
            if nargin == 5 % An input function u is set
                e(1) = e(1) + phi*this.M2*ut + ut'*this.M3*ut;    
            end
            this.neg_e1 = this.neg_e1 || e(1) < 0;
            
            e(1) = sqrt(abs(e(1)));
            %e(1) = sqrt(max(e(1),0));
            
            %% Normal computations
            if ~isempty(this.ReducedModel.FullModel.Approx)
                a = this.ReducedModel.FullModel.Approx;
            else
                a = this.ReducedModel.FullModel.System.f;
            end
            %idx = find(abs(this.xfull(1,:) - t) < eps,1);
            idx = this.cnt;
            xf = this.xfull(2:end,idx);
            fax = this.fax(:,idx);
            Vz = this.ReducedModel.V*x(1:end-this.ExtraODEDims);
            fVz = a.evaluate(Vz,t,mu);

            beta = 0;
            if this.Version == 1
                xdiff = sum((xf-Vz).^2);
                if xdiff ~= 0
                    beta = sqrt(sum((fax-fVz).^2)/xdiff); %#ok<*PROP>
                else
                    beta = 0;
                end
            elseif this.Version == 2
                if size(xf,1) < 300
                    Jx = a.getStateJacobian(xf,t,mu);
                    JVz = a.getStateJacobian(Vz,t,mu);
                    beta = max(norm(Jx),norm(JVz));
                end
            % Comparison lemma
            %e(2) = beta3*eold(2) + e(1);
            end
            e(2) = beta*(eold(1) + eold(2));
            this.betas(:,end+1) = [t; beta];
            this.cnt = this.cnt + 1;
        end
        
        function e0 = getE0(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = [this.ReducedModel.getExo(mu); 0];%; 0; 0; this.ReducedModel.getExo(mu)];
            
            %%%%%%%%%%%%% Experimental part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            rm = this.ReducedModel;
            [t, x] = rm.FullModel.computeTrajectory(mu, rm.System.inputidx);
            this.xfull = [t; x];
            if size(x,1) >= 300
                warning('a:b','Main dimension larger than 300. Disabling jacobian estimation version.');
            end
            if ~isempty(rm.FullModel.Approx)
                this.fax = rm.FullModel.Approx.evaluate(x,t,mu);
            else
                this.fax = rm.FullModel.System.f.evaluate(x,t,mu);
            end
        end
    end
    
    methods(Access=protected)
        function postprocess(this, t, x, mu, inputidx)%#ok
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('CompWiseErrorEstimator:process','Integral part is all zero. Attention!');
            end
            this.StateError = eint(1,:) + eint(2,:);
            
            if this.neg_e1
                disp('ExpensiveKernelEstimator: Negative alpha(t) norms occurred. Used zero instead.');
                this.neg_e1 = false;
            end
        end
    end
    
end