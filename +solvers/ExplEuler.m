classdef ExplEuler < solvers.BaseSolver
    % Explicit forward euler ODE solver
    %
    % This solver uses the MaxStep property as timestep to be most
    % efficient. Set small enough for precision.
    %
    % Method description:
    % `` x_{i+1} = x_i + \Delta t f(t_i,x_i)``
    %
    % See also: solvers BaseSolver Heun
    
    methods
        
        function this = ExplEuler(MaxStep)
            % Constructor for the explicit euler solver.
            %
            % Specifying the dt parameter is optional, though not
            % specifying will generate a warning that it should be set
            % manually for each case.
            %
            % Parameters:
            % MaxStep: Maximum time step. Optional.
            this.Name = 'Explicit forward euler';
            if nargin == 1
                this.MaxStep = MaxStep;
            end
%             if nargin == 0
%                 warning('solvers:ExplEuler:No_dt_given','Explicit solvers should get a time stepsize');
%                 this.MaxStep = 0.05;
%             else
%                 this.MaxStep = dt;
%             end
        end
        
        function [tout,y] = solve(this, odefun, t, x0)
            % Get computation times
            [times, tout, outputtimes] = this.getCompTimes(t);
            % Initialize vector
            steps = length(times);
            y = [x0 zeros(length(x0),steps-1)];
            dt = times(2:end)-times(1:end-1);
            % Solve for each time step
            for idx = 2:steps;
                y(:,idx) = y(:,idx-1) + dt(idx-1)*odefun(times(idx-1),y(:,idx-1));
            end
            % Extract wanted values
            y = y(:,outputtimes);
        end
    end
    
end

