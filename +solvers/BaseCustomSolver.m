classdef BaseCustomSolver < solvers.BaseSolver
    % BaseCustomSolver: Base class for all self-implemented solvers.
    %
    % Adds another abstraction layer to distinguish between self-implemented ode solvers and third-party
    % solvers.
    %
    % @author Daniel Wirtz @date 2011-05-31
    %
    % @change{0,7,ts,2013-01-24} accelerated method getCompTimes
    %
    % @new{0,4,dw,2011-05-31} Added this class.
    %
    % @todo checkout code at
    % http://www.mymathlib.com/diffeq/embedded_runge_kutta/embedded_prince_dormand_v1_4_5.html
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    methods
        function this = BaseCustomSolver
            this = this@solvers.BaseSolver;
            this.Name = 'Base Custom Solver';
        end
        
        function [t, x] = solve(this, odefun, t, x0)
            % Solves the ode and wraps the actual solving method.
            %
            % Parameters:
            % odefun: A function handle to the ode function, satisfying the
            % interface also required by matlab's explicit ode solvers.
            % @type function_handle
            % t: The desired times `t_0,\ldots,t_N` as row vector. @type rowvec
            % x0: The initial value `x(0) = x_0` for `t=0` @type colvec
            %
            % Return values:
            % t: The times `t_i` @type rowvec
            % x: The solution of the ode at the time steps `t_0,\ldots,t_N`
            % as matrix. @type matrix
            
            % Get computation times
            %[times, outputtimes] = this.getCompTimes(t);
            [times, outidx] = this.getCompTimes(t);
            
            % Fire PreSolve event
            ed = solvers.SolverEventData(times);
            notify(this,'PreSolve',ed);
            
            x = this.customSolve(odefun, times, x0, outidx);
            
            % Fire PostSolve event
            ed.States = x;
            notify(this,'PostSolve',ed);
            
            % Extract wanted values
            t = times(outidx);
            % New: x is supposed to be returned in the size of effectively desired timesteps
            % (memory issues for large, long simulations)
            %x = x(:,outputtimes);
        end
    end
    
    methods(Access=private)
        function [times, outidx] = getCompTimes(this, t)
            % Computes the computation and effective output times for a
            % given input time vector t.
            %
            % This method is thought of as a convenience method for custom solver implementations
            % that merges the desired evaluation times with the times resulting from a possible
            % MaxStep setting.
            %
            % As maximum execution speed is wanted, MaxStep is used as `dt`
            % time-step size. So in order to set `dt` use the MaxStep
            % property.
            % This method returns the given values if MaxStep is empty.
            %
            % If 't' was an array, we have @code t =
            % times(outputtimes)@endcode for the resulting parameters.
            % Otherwise, if 't' was a two-element vector (start- and
            % endtime) we get @code times = t(1):this.MaxStep:t(2) @endcode
            %
            % Parameters:
            % t: The desired times t. Either a two element vector
            % containing start time `t_0` and end time `T` or a row vector
            % `t=[t_0, t_1, \ldots, T]`. @type rowvec<double>
            %
            % Return values:
            % times: The actual times at which to compute the solution.
            % outidx: The indices of times the output is actually desired
            % at, i.e. times(outidx) = t. @type rowvec<integer> @default
            % 1:length(t)
            %
            % @todo InitialStep mit einbauen!
            %
            % @change{0,3,dw,2011-04-14} Fixed a problem that would cause useless long computation
            % times that occured when a MaxStep is set but is larger than any time-step passed in
            % the t vector.
            %
            % @change{0,2,dw,2011-03-11} Removed 'tout' from the return
            % parameters since it can be computed either way via 'tout =
            % times(outputtimes)'.
            
            % Validity checks
            if any(abs(sort(t) - t) > 100*eps)
                error('The time vector must be strictly monotoneously increasing.');
            end
            
            % Default values
            % outputtimes = true(1,length(t));
            outidx = 1:length(t);
            times = t;
            if ~isempty(this.MaxStep)
                if numel(t) == 2
                    times = t(1):this.MaxStep:t(2);
                    % outputtimes = true(1,length(times));
                    outidx = [1,length(times)];
                else
                    len = length(t);
                    n = max(round((t(2:len)-t(1:len-1))/this.MaxStep),1);  % number of steps in interval [t(i-1),t(i)] after refinemant. If t(i-1)-t(i)>= 1.5*MaxStep, interval will be refined
                    m = ones(1,len);      % m(i) is the index of the originally i-th entry in t in times 
                    m(2:len) = cumsum(n) + m(2:len);% - (0:len-2);
                    times = zeros(1,m(len));
                    % outputtimes = times;
                    for i = 1:len-1
                        times(m(i):m(i+1))=  t(i):(t(i+1)-t(i))/(n(i)):t(i+1); %linspace(t(i),t(i+1),n(i));
                        % outputtimes(m(i)) = 1;
                        outidx(i) = m(i);
                    end
                    %outputtimes(end) = 1;
                    outidx(end) = length(times);
                end
                %outputtimes = logical(outputtimes);
                
                %                     % Find refinement indices
                %                     idx = fliplr(find(t(2:end)-t(1:end-1)-this.MaxStep>100*eps));
                %                     % If any "gaps" are present, fill up with MaxStep's
                %
                %                     if ~isempty(idx)
                %                         for i=idx
                %                             add = times(i):this.MaxStep:times(i+1);
                %                             times = [times(1:i-1) add times(i+2:end)];
                %                             outputtimes = [outputtimes(1:i) false(1,length(add)-2) outputtimes(i+1:end)];
                %                         end
                %                         tout = times(outputtimes);
                %                         if numel(tout) ~= numel(t) || any(abs(tout - t) > 100*eps)
                %                             error('Unexpected error: Computed time vector differs from the desired one.');
                %                             t%#ok
                %                         end
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % The abstract solve function for custom ODE solver implementations.
        %
        % In difference to the solvers.BaseSolver.solve method the times passed here are already
        % matched to possible MaxStep restrictions (thus extended).
        %
        % Parameters:
        % odefun: A function handle to the ode's function. Signature is
        % 'odefun(t,x)'
        % t: Either a two dimensional vector with t(1) < t(2) specifiying
        % the start and end time, or a greater or equal to three
        % dimensional, strictly monotoneously increasing vector explicitly
        % setting the desired output times. Depending on the MaxStep
        % property, the solver can work with a finer time step internally.
        % x0: The initial value
        % outputtimes: A logical vector indicating the output times of
        % effective interest.
        %
        % See also: getCompTimes
        %
        % @change{0,7,dw,2013-01-11} Passing the outputtimes parameter to the customSolve
        % method in order to allow for more memory-efficient implementations.
        x = customSolve(this, odefun, t, x0, outputtimes);
    end
    
    events
        % Gets fired before the solver starts.
        %
        % Carries an SolverEventData instance with the Times field set to the actual times being
        % used during solving.
        %
        % @attention This event is not fired when using a Matlab builtin solver.
        %
        % See also: solver.ode.SolverEventData
        PreSolve;
        
        % Gets fired after the solver has finished.
        %
        % Carries an SolverEventData instance with the Times field set to the actual times being
        % used during solving and the States field set to the corresponding state variable values.
        %
        % @attention This event is not fired when using a Matlab builtin solver.
        %
        % See also: solver.ode.SolverEventData
        PostSolve;
    end
    
end