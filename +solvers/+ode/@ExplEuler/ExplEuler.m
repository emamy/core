classdef ExplEuler < solvers.ode.BaseCustomSolver
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
    % See also: solvers BaseSolver BaseCustomSolver Heun
    %
    % @change{0,5,dw,2011-09-29} Added step-wise event implementation for real time
    % plotting.
    %
    % @change{0,4,dw,2011-05-31} Added a new middle class solvers.ode.BaseCustomSolver which
    % extracts the getCompTimes into a new abstraction layer.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @new{0,2,dw,2011-03-11} Added a c/mex implementation of the
    % algorithm. Turns out it is double the times slower than the matlab
    % native code version, so leaving it only in there for speed test
    % purposes (solvers.ode.ExplEuler.test_solveMex).
    
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
            this = this@solvers.ode.BaseCustomSolver;
            
            this.Name = 'Explicit forward euler';
            if nargin == 1
                this.MaxStep = MaxStep;
            end
        end
    end
    
    methods(Access=protected,Sealed)
        function x = customSolve(this, odefun, t, x0)%#ok
            % Solves the ODE
            %
            % There is a mex-implementation of the explicit euler scheme,
            % unfortunately it is almost double the times slower than
            % matlab! (So mexing files to get rid of for loops does not
            % work well :-))
            
            % Initialize vector
            steps = length(t);
            x = [x0 zeros(length(x0),steps-1)];
            dt = t(2:end)-t(1:end-1);
            
            % Solve for each time step
            ed = solvers.ode.SolverEventData;
            for idx = 2:steps;
                x(:,idx) = x(:,idx-1) + dt(idx-1)*odefun(t(idx-1),x(:,idx-1));
                ed.Times = t(idx);
                ed.States = x(:,idx);
                this.notify('StepPerformed',ed);
            end

            %y = this.solveMex(odefun, times, x0);
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
            s = solvers.ode.ExplEuler;
            
            times = 0:.05:1;
            iter = 100;
            tmex = zeros(1,iter);
            tmat = zeros(1,iter);
            
            for it = 1:iter
            
            t = tic;
            s.solveMex(@odefun, times, rand(10,1));
            tmex(it) = toc(t);
            
            t = tic;
            s.solve(@odefun, times, rand(10,1));
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

