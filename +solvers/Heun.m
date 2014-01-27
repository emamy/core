classdef Heun < solvers.BaseCustomSolver
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
% @change{0,5,dw,2011-10-16} Adopted to the new BaseSolver.RealTimeMode flag.
%
% @change{0,5,dw,2011-09-29} Added step-wise event implementation for real time
% plotting.
%
% @change{0,4,dw,2011-05-31} Added a new middle class solvers.BaseCustomSolver which
% extracts the getCompTimes into a new abstraction layer.
%
% @todo cater for time-constant mass matrices (get LU once before)
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
        
    methods
        
        function this = Heun(MaxStep)
            % Constructor for the explicit euler solver.
            %
            % Specifying the dt parameter is optional, though not
            % specifying will generate a warning that it should be set
            % manually for each case.
            %
            % Parameters:
            % MaxStep: Maximum time step. @default [] @type double
            
            this = this@solvers.BaseCustomSolver;
            
            this.Name = 'Explicit Heun''s method';
            this.SolverType = solvers.SolverTypes.Explicit;
            if nargin == 1
                this.MaxStep = MaxStep;
            end
        end
    end
    
    methods(Access=protected,Sealed)
        function x = customSolve(this, odefun, t, x0, outputtimes)
            % Solves the ode using Heuns method.
            %
            % Parameters:
            % odefun: A function handle to the ode function, satisfying the
            % interface also required by matlab's explicit ode solvers.
            % @type function_handle
            % t: The desired times `t_0,\ldots,t_N` as row vector. @type rowvec
            % x0: The initial value `x(0) = x_0` for `t=0` @type colvec
            % outputtimes: index vector indicating at which of the times in t
            % the output is actually desired. The solution will be returned
            % at `t_0,\ldots,t_N` = t(outputtimes).
            % @type rowvec<integer>
            %
            % Return values:
            % x: The solution of the ode at the time steps `t_0,\ldots,t_N`
            % as matrix. @type matrix
            
            % Initialize result
            steps = length(t);
            dt = t(2:end)-t(1:end-1);
            
            rtm = this.RealTimeMode;
%             if rtm
%                 ed = solvers.SolverEventData;
%                 x = [];
%             else
%                 x = [x0 zeros(size(x0,1),steps-1)];
%             end
            % Initialize output index counter
            outidx = 2;
            if rtm
                ed = solvers.SolverEventData;
                x = [];
            else
                effsteps = length(outputtimes);
                % Create return matrix in size of effectively desired timesteps
                x = [x0 zeros(size(x0,1),effsteps-1)];
            end
            % Solve for each time step
            oldx = x0;
            for idx = 2:steps;
                f = odefun(t(idx-1),oldx);
                hlp = (dt(idx-1)/2)*(f + odefun(t(idx),oldx + dt(idx-1)*f));
                
                % Check if a mass matrix is present
                if ~isempty(this.M)
%                     newx = oldx + this.M.evaluate(t(idx-1))\hlp;
                    [L,U] = this.M.getLU(t(idx-1)); 
                    newx = U\(L\(this.M.evaluate(t(idx-1))*oldx + hlp));
                else
                    newx = oldx + hlp;
                end
                
%                 if rtm
%                     ed.Times = t(idx);
%                     ed.States = newx;
%                     this.notify('StepPerformed',ed);
%                 else
%                     x(:,idx) = newx;%#ok
%                 end
                % Only produce output at wanted timesteps
                if outputtimes(outidx) == idx
                    if rtm                        
                        % Real time mode: Fire StepPerformed event
                        ed.Times = t(idx);
                        ed.States = newx;
                        this.notify('StepPerformed',ed);
                        % Normal mode: Collect solution in result vector
                    else
                        x(:,outidx) = newx;
                    end
                    outidx = outidx+1;
                end
                oldx = newx;
            end
        end
    end
    
end

