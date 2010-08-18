classdef BaseSolver < handle
    % Base class for all KerMor ODE solvers
    %
    % Simply defines an interfaces for a solve function and provides common
    % ODE solver properties.
    %
    % @todo: Write tests for solvers.
    
    properties
        % Maximum time step for solver.
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
        function [times, tout, outputtimes] = getCompTimes(this, t, dt)
            % Computes the computation and effective output times for a
            % given input time vector t and the desired timestep dt.
            %
            % Parameters:
            % t: The desired times t. Either a two element vector
            % containing start time `t_0` and end time `T` or a row vector
            % `t=[t_0, t_1, \ldots, T]`.
            % dt: The internal timestep to use during computation.
            %
            % Return values:
            % times: The actual times at which to compute the solution.
            % tout: The effective times that are being returned along with
            % the solution values.
            % outputtimes: A logical row vector of the size of times that
            % indicates which element from the times vector also belongs to
            % the tout vector.
            
            if numel(t) == 2
                times = t(1):dt:t(2);
                outputtimes = true(1,length(times));
                tout = times;
            else
                % Else: Mix together full detailed times with
                full = t(1):dt:t(end);
                times = union(t,full);
                outputtimes = false(1,length(times));
                [xx,tidx,xxx] = intersect(times,t);%#ok
                outputtimes(tidx) = true;
                tout = t;
            end
        end
    end
    
    methods(Abstract)
        % The abstract solve function for the ODE solver.
        [t,y] = solve(odefun, t, x0, opts);
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
        
        function res = testSolverSpeedTest
            m = models.synth.KernelTest1(1400);
            
            perform(solvers.MLWrapper(@ode23));
            perform(solvers.MLWrapper(@ode45));
            perform(solvers.ExplEuler);
            perform(solvers.Heun);
            
            res = 1;
            
            function perform(solver)
                m.ODESolver = solver;
                tic;
                m.offlineGenerations;
                r = m.buildReducedModel;
                t = toc;
                fprintf('Using solver %s\n',m.ODESolver.Name);
                fprintf('Offline generations time: %f\n',t);
                [ti,x,xr,t,tr,tr_noerr] = r.getTrajectories;
                fprintf('Online simulations time\nFull detail: %fs\nReduced with error estimator: %fs\nReduced without error estimation:%fs\n\n',t,tr,tr_noerr);
            end
        end
        
    end
    
end

