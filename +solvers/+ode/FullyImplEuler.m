classdef FullyImplEuler < solvers.ode.BaseCustomSolver & solvers.ode.AImplSolver
% FullyImplSolver: Solver for fully nonlinear ODE's (using Newton iterations)
%
%
%
% @author Daniel Wirtz @date 2012-02-03
%
% @new{0,6,dw,2012-02-03} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(Access=private)
        % The model which the solver is used in
        model;
    end
    
    methods
        function this = FullyImplEuler(model)
            if ~isa(model.System.f,'dscomponents.IJacobian')
                error('FullyImplEuler can only be used with system functions that implement dscomponents.IJacobian');
            end
            this.model = model;
            % "Disable" MaxStep DPCM warning as implicit solvers are stable
            this.MaxStep = [];
        end
    end
    
    methods(Access=protected)
        function x = customSolve(this, odefun, t, x0)
           
            % Checks
            steps = length(t);
            dt = t(2)-t(1);
            if any(t(2:end)-t(1:end-1) - dt > 100*eps)
                error('non-equidistant dt timesteps.');
            end
            
            % Initializations
            rtm = this.RealTimeMode;
            if rtm
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                x = [x0 zeros(size(x0,1),steps-1)];
            end
            s = this.model.System;
            
            % Check if a mass matrix is present, otherwise assume identity matrix
            M = speye(length(x0));
            mdep = false;
            if ~isempty(s.M)
                mdep = s.M.TimeDependent;
                if ~mdep
                    M = s.M.evaluate(0); 
                end
            end
                                    
            % Solve for each time step
            oldx = x0;
            for idx = 1:steps-1;
                % Case: time-dependent Mass Matrix
                if mdep
                    % Question: choose t(idx) or t(idx+1)?
                    % Answer (Chris): t(idx+1), because this choice is independent of time-discr. of d/dt x !
                    % and for impl Euler the evaluation time is the "future"-timestep
                    M = s.M.evaluate(t(idx+1));
                end
                
                %% Implementation part:
                % principal equation: M(t)x'(t) = g(x(t),t)
                % discretized M(t)x(t+) = M(t)x(t) + dt*g(x(t+),t+), t+ = t+\delta t
                
                % g is given by odefun handle: g(x(t),t) = odefun(<t>,<x>)
                % \delta t is given by "dt" (above)
                % \Nabla g is accessible via s.f.getStateJacobian(<x>,<t>,s.mu)
                
                % write result of newton iteration to "newx"
                
                %% TODO: Newton-Iteration 'by hand'
%                 for newton_it = 1:max_newton_steps
%                 
%                     delta_x = (M - dt * s.f.getStateJacobian(oldx, t(idx+1), s.mu)) \ ( M*x )
%                     
%                 end
                %% Matlab fsolve
                dis = 'iter';
                %dis = 'final-detailed';
                options_fsolve = optimset( 'Display', dis, 'Jacobian',...
                        'on', 'MaxIter', 2000, 'MaxFunEvals', 10000, 'TolFun', 1e-3);
                nonlin_fun = @(x) deal(M * (x - oldx) - dt*odefun(t(idx+1), x),...
                                   M - dt * s.f.getStateJacobian(x, t(idx+1), s.mu));
                % Nullstelle der schwachen Form finden
                [newx, fval, exitflag, output, J_check] = fsolve( nonlin_fun, oldx, options_fsolve );
                
                %% Postprocessing
                % Real time mode: Fire StepPerformed event
                if rtm
                    ed.Times = t(idx);
                    ed.States = newx;
                    this.notify('StepPerformed',ed);
                % Normal mode: Collect solution in result vector
                else
                    x(:,idx) = newx;%#ok
                end
                oldx = newx;
            end
        end
    end
    
end