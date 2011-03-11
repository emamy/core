classdef BaseSolver < handle
    % Base class for all KerMor ODE solvers
    %
    % Simply defines an interfaces for a solve function and provides common
    % ODE solver properties.
    %
    % @todo Write tests for solvers.
    % @todo Rename BaseSolver to BaseODESolver
    
    properties
        % Maximum time step for solver.
        %
        % Upon time-step computation via getCompTimes this value is used as
        % `dt` inner time-step.
        %
        % Default: [] (Use solver's default settings)
        MaxStep = [];
        
        % The initial step size for the solver.
        % So far only used by the MLWrapper class for the native matlab
        % solvers.
        %
        % Default: [] (User solver's default settings)
        InitialStep = [];
    end
    
    properties(SetAccess=protected)
        % The solver's name
        Name = 'KerMor BaseSolver';
    end
    
    methods(Access=protected)
        function [times, outputtimes] = getCompTimes(this, t)
            % Computes the computation and effective output times for a
            % given input time vector t. 
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
            % `t=[t_0, t_1, \ldots, T]`.
            %
            % Return values:
            % times: The actual times at which to compute the solution.
            % outputtimes: A logical row vector of the size of times that
            % indicates which element from the times vector also belongs to
            % the effective output times vector.
            %
            % @todo InitialStep mit einbauen!
            %
            % @change{0,3,dw,2011-03-11} Removed 'tout' from the return
            % parameters since it can be computed either way via 'tout =
            % times(outputtimes)'.
            
            % Validity checks
            if any(abs(sort(t) - t) > 100*eps)
                error('The time vector must be strictly monotoneously increasing.');
            end
            
            % Default values
            outputtimes = true(1,length(t));
            times = t;
            if ~isempty(this.MaxStep)
                if numel(t) == 2
                    times = t(1):this.MaxStep:t(2);
                    outputtimes = true(1,length(times));
                else    
                    % Find refinement indices
                    idx = fliplr(find(abs(t(2:end)-t(1:end-1)-this.MaxStep)>100*eps));
                    % If any "gaps" are present, fill up with MaxStep's
                    if ~isempty(idx)
                        for i=idx
                            add = times(i):this.MaxStep:times(i+1);
                            times = [times(1:i-1) add times(i+1:end)];
                            outputtimes = [outputtimes(1:i) false(1,length(add)-1) outputtimes(i+1:end)];
                        end
                        tout = times(outputtimes);
                        if numel(tout) ~= numel(t) || any(abs(tout - t) > 100*eps)
                            error('Unexpected error: Computed time vector differs from the desired one.');
                            t%#ok
                        end
                    end
                end
            end
        end
    end
    
    methods(Abstract)
        % The abstract solve function for the ODE solver.
        %
        % Parameters:
        % odefun: A function handle to the ode's function. Signature is
        % 'odefun(t,x)'
        % t: Either a two dimensional vector with t(1) < t(2) specifiying
        % the start and end time, or a greater or equal to three
        % dimensional, strictly monotoneously increasing vector explicitly
        % setting the desired output times. Depending on the MaxStep
        % property, the solver can work with a finer time step internally.
        solve(this, odefun, t, x0, opts);
        
        %[t,y] = solve(odefun, t, x0, opts);
    end
    
    methods
        function set.MaxStep(this, value)
            if isempty(value)
                this.MaxStep = value;
                return;
            elseif ~isposrealscalar(value)
                error('Positive real scalar expected.');
            elseif value == 0
                error('Maximum time step must be greater than zero. Use [] to unset.');
            end
            this.MaxStep = value;
        end
        
        function set.InitialStep(this, value)
            if ~isposrealscalar(value)
                error('Positive real scalar expected.');
            end
            this.InitialStep = value;
        end
    end
    
    methods(Static)
        
        function res = test_SolverSpeedTest
            m = models.synth.KernelTest(200);
            perform(solvers.ExplEuler);
            perform(solvers.MLWrapper(@ode23));
            perform(solvers.MLWrapper(@ode45));
            perform(solvers.Heun);
            
            res = 1;
            
            function perform(solver)
                m.ODESolver = solver;
                tic;
                m.offlineGenerations;
                r = m.buildReducedModel;
                r.ErrorEstimator.Iterations = 0;
                t = toc;
                fprintf('Using solver %s\n',m.ODESolver.Name);
                fprintf('Offline generations time: %f\n',t);
                [ti,x,xr,t,tr,tr_noerr] = r.getTrajectories;
                fprintf('Online simulations time\nFull detail: %fs\nReduced with error estimator: %fs\nReduced without error estimation:%fs\n\n',t,tr,tr_noerr);
            end
        end
        
    end
    
end

