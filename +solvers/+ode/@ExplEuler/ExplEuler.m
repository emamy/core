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
% @change{0,5,dw,2011-10-16} Adopted to the new BaseSolver.RealTimeMode flag.
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
%
% See also: solvers BaseSolver BaseCustomSolver Heun
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        
        function this = ExplEuler(MaxStep)
            % Constructor for the explicit euler solver.
            %
            % Specifying the dt parameter is optional, though not
            % specifying will generate a warning that it should be set
            % manually for each case.
            %
            % Parameters:
            % MaxStep: Maximum time step. @default [] @type double
            this = this@solvers.ode.BaseCustomSolver;
            
            this.Name = 'Explicit forward euler';
            if nargin == 1
                this.MaxStep = MaxStep;
            end
        end
    end
    
    methods(Access=protected,Sealed)
        function x = customSolve(this, odefun, t, x0)
            % Solves the ODE using the explicit Euler method.
            %
            % There is a mex-implementation of the explicit euler scheme,
            % unfortunately it is almost double the times slower than
            % matlab! (So mexing files to get rid of for loops does not
            % work well :-))
            %
            % Parameters:
            % odefun: A function handle to the ode function, satisfying the
            % interface also required by matlab's explicit ode solvers.
            % @type function_handle
            % t: The desired times `t_0,\ldots,t_N` as row vector. @type rowvec
            % x0: The initial value `x(0) = x_0` for `t=0` @type colvec
            %
            % Return values:
            % x: The solution of the ode at the time steps `t_0,\ldots,t_N`
            % as matrix. @type matrix
            
            % Initialize result
            steps = length(t);
            dt = t(2:end)-t(1:end-1);
            
            rtm = this.RealTimeMode;
            if rtm
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                x = [x0 zeros(size(x0,1),steps-1)];
            end
            
            % Solve for each time step
            oldx = x0;
            for idx = 2:steps;
                
                hlp = dt(idx-1)*odefun(t(idx-1),oldx);
                % Check if a mass matrix is present
                if ~isempty(this.M)
                    newx = oldx + this.M.evaluate(t(idx-1))\hlp;
%                     [L,U,Q,P] = this.M.getLU(t(idx-1)); 
%                     newx = oldx + Q*(U\(L\(P*hlp)));
%                     newx = oldx + U\(L\hlp);
                else
                    newx = oldx + hlp;
                end
                
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

