classdef ReducedSecondOrderSystem < models.ReducedSystem & models.BaseSecondOrderSystem
    %REDUCEDSECONDORDERSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function this = ReducedSecondOrderSystem(model)
            % Creates a new base dynamical system class instance.
            this = this@models.BaseSecondOrderSystem(model);
        end
        
        function setReducedModel(this, rmodel)
            setReducedModel@models.ReducedSystem(this, rmodel);
            % Additional steps
        end
        
        function odefun = getODEFun(this)
            % Determine correct ODE function (A,f,B combination)
            xarg = 'x';
            est = this.Model.ErrorEstimator;
            haveest = ~isempty(est) && est.Enabled;
            if haveest
                xarg = 'x(1:end-est.ExtraODEDims,:)';
            end
            
            str = {};
            if ~isempty(this.A)
                str{end+1} = sprintf('this.A.evaluate(%s, t)',xarg);
            end
            if ~isempty(this.f)
                str{end+1} = sprintf('this.f.evaluate(%s, t)',xarg);
            end
            if ~isempty(this.B) && ~isempty(this.inputidx)
                str{end+1} = 'this.B.evaluate(t, this.mu)*this.u(t)';
            end
            funstr = Utils.implode(str,' + ');
            if haveest
                funstr = ['[' funstr '; est.evalODEPart(x, t, this.u(t))]'];
            end
            odefun = eval(['@(t,x)' funstr]);
        end
        
        function x0 = getX0(this, mu)
            % Gets the initial value of `x_0(\mu)`.
            %
            % If the estimator is enabled, x0 is extended by the e0
            % components of the error estimator.
            
            x0 = getX0@models.BaseFirstOrderSystem(this, mu);
            m = this.Model;
            if ~isempty(m.ErrorEstimator) && m.ErrorEstimator.Enabled
                x0 = [x0; m.ErrorEstimator.getE0(mu)];
            end
        end
        
    end
    
    methods(Access=protected)
        function val = getDerivativeDirichletValues(this, t)
            val = getDerivativeDirichletValues@models.BaseSecondOrderSystem(this, t);
        end
    end
    
end

