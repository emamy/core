classdef LinearImplEuler < solvers.ode.BaseCustomSolver & solvers.ode.AImplSolver
% LinearImplEuler: Implicit euler solver for linear core functions
%
% This class is a "hack" since it explicitly accesses the model class and
% uses the system's components `f`, `B` and `u(t)` directly instead of
% using the provided odefun handle. Due to this, it is not possible to use
% this solver in conjunction with any error estimators as long as they
% extend the reduced system by any dimension.
%
% @author Daniel Wirtz @date 2011-12-06
%
% @new{0,6,dw,2011-12-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        model;
    end
    
    methods
        
        function this = LinearImplEuler(model)
            this.model = model;
            % "Disable" MaxStep DPCM warning as implicit solvers are stable
            this.MaxStep = [];
        end
    end
    
    methods(Access=protected)
        function x = customSolve(this, ~, t, x0)
            s = this.model.System;
            if ~isa(s.f,'dscomponents.LinearCoreFun') && ~isa(s.f,'dscomponents.AffLinCoreFun')
                error('LinearImplEuler only works for dscomponents.{LinearCoreFun|AffLinCoreFun} (subclasses).');
            end
            if isa(s.Model,'models.ReducedModel') && ~isempty(s.Model.ErrorEstimator)...
                && s.Model.ErrorEstimator.Enabled && s.Model.ErrorEstimator.ExtraODEDims > 0
                error('Cannot use this solver with reduced models that have an error estimator enabled with ExtraODEDims > 0 (this solver overrides the odefun handle)');
            end
            % Initialize result
            steps = length(t);
            dt = t(2)-t(1);
            if any(t(2:end)-t(1:end-1) - dt > 100*eps)
                error('non-equidistant dt timesteps.');
            end
            
            rtm = this.RealTimeMode;
            if rtm
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                x = [x0 zeros(size(x0,1),steps-1)];
            end
            
            % Check if a mass matrix is present
            if ~isempty(s.M)
                M = s.M.evaluate(0); 
            else
                M = eye(size(s.f.A));
            end
            if isa(s.f,'dscomponents.AffLinCoreFun')
                % Warning: no time-dependence implemented here
                A = s.f.AffParamMatrix.compose(0, s.mu);
            else
                A = s.f.A;
            end
            %[l,u] = lu(M + dt * A);
            Ai = inv(M + dt * A);
            
            % Solve for each time step
            oldx = x0;
            for idx = 2:steps;
                RHS = M*oldx + dt*s.f.b;
                if ~isempty(s.u)
                    RHS = RHS + dt*s.B.evaluate(t, s.mu)*s.u(t);
                end
                %newx = u\(l\RHS);
                newx = Ai * RHS;
                %newx = (M + dt * A)\RHS;
                
                if rtm
                    ed.Times = t(idx);
                    ed.States = newx;
                    this.notify('StepPerformed',ed);
                else
                    x(:,idx) = newx;%#ok
                end
                oldx = newx;
            end
        end
    end
    
end