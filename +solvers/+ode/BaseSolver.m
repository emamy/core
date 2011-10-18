classdef BaseSolver < KerMorObject
    % Base class for all KerMor ODE solvers
    %
    % Simply defines an interfaces for a solve function and provides common
    % ODE solver properties.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @todo Write tests for solvers.
    
    properties(SetObservable)
        % Maximum time step for solver.
        %
        % Upon time-step computation via getCompTimes this value is used as
        % `dt` inner time-step.
        %
        % @propclass{critical} Too large steps due to high time-step distances in the passed times
        % vector `t` may lead to errorneous results. This property limits the maximum time-step size
        % used in the implementations. Set to [] in order to rely on the times `t`.
        %
        % @default [] (Use solver's default settings)
        MaxStep = [];
        
        % The initial step size for the solver.
        % So far only used by the MLWrapper class for the native matlab
        % solvers.
        %
        % @propclass{important} Some odes require a certain (small) initial step. 
        %
        % @default [] (User solver's default settings)
        InitialStep = [];
    end
    
    properties(SetAccess=protected)
        % The solver's name
        %
        % @propclass{data} Define your own solver name!
        Name = 'KerMor BaseSolver class';
    end
    
    methods
        function this = BaseSolver
            this = this@KerMorObject;
            this.registerProps('MaxStep','InitialStep','Name');
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
        % @type rowvec
        % x0: The initial value @type colvec
        %
        % Return values:
        % t: The times `t_i` which equal the given times or the linspace
        % between both given times (two-elem t param) with distance `dt`
        % @type rowvec
        % x: The solution of the ode at times `t_i` @type matrix
        %
        [t,x] = solve(this, odefun, t, x0);
    end
    
    %% Getter & Setter
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
    
    events
        % Gets fired when an ODE solver performs an intermediate step
        %
        % @see models.BaseModel.RealTimePlotting
        StepPerformed;
    end
    
    methods(Static)
        
        function res = test_SolverSpeedTest
            m = models.synth.KernelTest(200);
            perform(solvers.ode.ExplEuler);
            perform(solvers.ode.MLWrapper(@ode23));
            perform(solvers.ode.MLWrapper(@ode45));
            perform(solvers.ode.Heun);
            
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

