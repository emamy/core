classdef ExpensiveBetaEstimator < error.BaseKernelEstimator
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
        xfull;
        fax;
        cnt;
    end
    
    methods
        function this = ExpensiveBetaEstimator(rmodel)
            this = this@error.BaseKernelEstimator;
            if nargin == 1
                this.setReducedModel(rmodel);
            end
            this.cnt = 1;
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.ExpensiveBetaEstimator;
            copy = clone@error.BaseKernelEstimator(this, copy);
        end
        
        function clear(this)
            clear@error.BaseKernelEstimator(this);
            this.cnt = 1;
        end
        
        function e0 = init(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = init@error.BaseKernelEstimator(this, mu);
            
            %% Expensive computations
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
        
        function b = getBeta(this, x, t, mu)
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
            
            b = 0;
            if this.Version == 1
                xdiff = sum((xf-Vz).^2);
                if xdiff ~= 0
                    b = sqrt(sum((fax-fVz).^2)/xdiff); %#ok<*PROP>
                else
                    b = 0;
                end
            elseif this.Version == 2
                if size(xf,1) < 300
                    Jx = a.getStateJacobian(xf,t,mu);
                    JVz = a.getStateJacobian(Vz,t,mu);
                    b = max(norm(Jx),norm(JVz));
                end
            end
            this.cnt = this.cnt + 1;
        end
        
        function postprocess(this, t, x, mu, inputidx)%#ok
            this.StateError = x(end,:);
        end
        
    end
    
end