classdef Heun < solvers.ode.BaseSolver
    % ODE solver implementing the method of heun
    %
    % Method description:
    % `` x_{i+1} = x_i + \frac{\Delta t}{2}f(t_i,x_i) + \frac{\Delta
    % t}{2}f(t_{i+1},\nu)``
    % `` \nu = x_i + \Delta t f(t_i,x_i) ``
    %
    % See also: solvers BaseSolver ExplEuler
        
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
            
            this = this@solvers.ode.BaseSolver;
            
            this.Name = 'Explicit Heun''s method';
            if nargin == 1
                this.MaxStep = MaxStep;
            end
        end
        
        function [tout,y] = solve(this, odefun, t, x0)
            % Get computation times
            [times, outputtimes] = this.getCompTimes(t);
            % Initialize result
            steps = length(times);
            y = [x0 zeros(size(x0,1),steps-1)];
            dt = times(2:end)-times(1:end-1);
            % Solve for each time step
            for idx = 2:steps;
                f = odefun(times(idx-1),y(:,idx-1));
                hlp = y(:,idx-1) + dt(idx-1)*f;
                y(:,idx) = y(:,idx-1) + (dt(idx-1)/2)*(f + odefun(times(idx),hlp));
            end
            y = y(:,outputtimes);
            tout = times(outputtimes);
            
%             cur = x0; outidx = 2;
%             steps = length(times);
%             % Solve for each time step
%             for idx = 2:steps;
%                 f = odefun(times(idx-1),cur);
%                 hlp = cur + this.MaxStep*f;
%                 cur = cur + (this.MaxStep/2)*(f + odefun(times(idx),hlp));
%                 if outputtimes(idx)
%                     y(:,outidx) = cur;
%                     outidx = outidx+1;
%                 end
%             end            
        end
    end
    
end

