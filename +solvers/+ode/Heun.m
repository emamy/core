classdef Heun < solvers.ode.BaseCustomSolver
    % ODE solver implementing the method of heun
    %
    % Method description:
    % `` x_{i+1} = x_i + \frac{\Delta t}{2}f(t_i,x_i) + \frac{\Delta
    % t}{2}f(t_{i+1},\nu)``
    % `` \nu = x_i + \Delta t f(t_i,x_i) ``
    %
    % See also: solvers BaseSolver ExplEuler
    %
    % @author Daniel Wirtz @date 2010-11-03
    %
    % @change{0,4,dw,2011-05-31} Added a new middle class solvers.ode.BaseCustomSolver which
    % extracts the getCompTimes into a new abstraction layer.
        
    methods
        
        function this = Heun(MaxStep)
            % Constructor for the explicit euler solver.
            %
            % Specifying the dt parameter is optional, though not
            % specifying will generate a warning that it should be set
            % manually for each case.
            %
            % Parameters:
            % MaxStep: Maximum time step. Optional.
            
            this = this@solvers.ode.BaseCustomSolver;
            
            this.Name = 'Explicit Heun''s method';
            if nargin == 1
                this.MaxStep = MaxStep;
            end
        end
    end
    
    methods(Access=protected,Sealed)
        function x = customSolve(this, odefun, t, x0)%#ok
            % Initialize result
            steps = length(t);
            x = [x0 zeros(size(x0,1),steps-1)];
            dt = t(2:end)-t(1:end-1);
            % Solve for each time step
            for idx = 2:steps;
                f = odefun(t(idx-1),x(:,idx-1));
                hlp = x(:,idx-1) + dt(idx-1)*f;
                x(:,idx) = x(:,idx-1) + (dt(idx-1)/2)*(f + odefun(t(idx),hlp));
            end
        end
    end
    
end

