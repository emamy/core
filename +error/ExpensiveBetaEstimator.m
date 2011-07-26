classdef ExpensiveBetaEstimator < error.BaseCompLemmaEstimator
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
    
    properties(Access=private,Transient)
        xfull;
        fax;
        cnt;
    end
    
    methods
        function this = ExpensiveBetaEstimator(rmodel)
            this = this@error.BaseCompLemmaEstimator;
            if nargin == 1
                this.setReducedModel(rmodel);
            end
            this.cnt = 1;
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.ExpensiveBetaEstimator;
            copy = clone@error.BaseCompLemmaEstimator(this, copy);
        end
        
        function clear(this)
            clear@error.BaseCompLemmaEstimator(this);
            this.cnt = 1;
        end
        
        function prepareConstants(this, mu, inputidx)
            prepareConstants@error.BaseCompLemmaEstimator(this, mu, inputidx);
            
            rm = this.ReducedModel;
            [t, x] = rm.FullModel.computeTrajectory(mu, inputidx);
            this.xfull = [t; x];
            mu = repmat(mu,1,length(t));
            if ~isempty(rm.FullModel.Approx)
                this.fax = rm.FullModel.Approx.evaluate(x,t,mu);
            else
                this.fax = rm.FullModel.System.f.evaluate(x,t,mu);
            end
        end
    end
    
    methods(Access=protected)
        
        function b = getBeta(this, x, t, mu)
            % Throw away the last part, not needed here
            x = x(1:end-1);
            if ~isempty(this.ReducedModel.FullModel.Approx)
                a = this.ReducedModel.FullModel.Approx;
            else
                a = this.ReducedModel.FullModel.System.f;
            end
            %idx = find(abs(this.xfull(1,:) - t) < eps,1);
            idx = this.cnt;
            xf = this.xfull(2:end,idx);
            fax = this.fax(:,idx);
            if size(this.ReducedModel.V,1) ~= size(x,1)
                x = this.ReducedModel.V*x;
            end
            fx = a.evaluate(x,t,mu);
            
%             b = 0;
%             if this.Version == 1
                xdiff = sum((xf-x).^2);
                if xdiff ~= 0
                    b = sqrt(sum((fax-fx).^2)/xdiff); %#ok<*PROP>
                else
                    b = 0;
                end
%             elseif this.Version == 2
%                 if size(xf,1) < 300
%                     Jx = a.getStateJacobian(xf,t,mu);
%                     JVz = a.getStateJacobian(Vz,t,mu);
%                     b = max(norm(Jx),norm(JVz));
%                 end
%             end
            this.cnt = this.cnt + 1;
        end
        
        function postprocess(this, t, x, mu, inputidx)%#ok
            this.StateError = x(end,:);
        end
        
    end
    
end