classdef ExplEuler < solvers.BaseSolver
    % Explicit forward euler ODE solver
    %
    % This solver uses the MaxStep property as timestep to be most
    % efficient. Set small enough for precision.
    %
    % Method description:
    % `` x_{i+1} = x_i + \Delta t f(t_i,x_i)``
    %
    % @author Daniel Wirtz @date 12.03.2011
    %
    % See also: solvers BaseSolver Heun
    %
    % @new{0,2,dw,2011-03-11} Added a c/mex implementation of the
    % algorithm. Turns out it is double the times slower than the matlab
    % native code version, so leaving it only in there for speed test
    % purposes (solvers.ExplEuler.test_solveMex).
    
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
            % Solves the ODE
            %
            % There is a mex-implementation of the explicit euler scheme,
            % unfortunately it is almost double the times slower than
            % matlab! (So mexing files to get rid of for loops does not
            % work well :-))
            
            % Get computation times
            [times, outputtimes] = this.getCompTimes(t);
            
            % Initialize vector
            steps = length(times);
            y = [x0 zeros(length(x0),steps-1)];
            dt = times(2:end)-times(1:end-1);
            % Solve for each time step
            for idx = 2:steps;
                y(:,idx) = y(:,idx-1) + dt(idx-1)*odefun(times(idx-1),y(:,idx-1));
            end
            %y = this.solveMex(odefun, times, x0);
            
            % Extract wanted values
            y = y(:,outputtimes);
            tout = times(outputtimes);
        end
        
        
    end
    
    methods(Access=public)
        % brief upper
        %
        % detail
        [tout, y] = solveMex(this, odefun, t, x0);
        % brief lower
        %
        % detail lower
    end
    
    methods(Static)
        function res = test_solveMex
            s = solvers.ExplEuler;
            
            times = 0:.05:1;
            iter = 100;
            tmex = zeros(1,iter);
            tmat = zeros(1,iter);
            
            for it = 1:iter
            
            t = tic;
            [ti,y] = s.solveMex(@odefun, times, rand(10,1));
            tmex(it) = toc(t);
            
            t = tic;
            [ti,y] = s.solve(@odefun, times, rand(10,1));
            tmat(it) = toc(t);
            
            end
            
            fprintf('Mex time: %f, Mat time: %f\n',mean(tmex),mean(tmat));
            
            res = true;
            
            function y = odefun(t,x)
                y = x+t*x;
            end
        end
    end
end

