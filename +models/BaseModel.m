classdef BaseModel < KerMorObject
% BaseModel: Base class for both full and reduced models.
%
% This class gathers all common functionalities of models in the
% KerMor framework.
% The most important method would be @code [t,y] =
% simulate(mu,inputidx) @endcode which computes the system's solution
% for given `\mu` and input number (if applicable).  Also a plot
% wrapper is provided that refers to the plotting methods within the
% model's system.
%
% @author Daniel Wirtz @date 19.03.2010
%
% @change{0,6,dw,2012-05-26} Adopted way of computing the jacobian and
% jsparsity pattern for ode solvers according to new possibility of having
% a (affine) linear models.BaseDynSystem.A component.
%
% @change{0,6,dw,2011-12-14} Introduced a \c ctime argument to the simulate and
% computeTrajectory methods. This motivates from several caching strategies that may be applied
% which lead to less computation time due to simple trajectory lookup. however, when comparing
% error estimators this effect irritates the results since subsequent calls to trajectory
% computations in the context of different estimators lead to different computation times.
%
% @change{0,6,dw,2011-11-25} Made T a dependent property and added consistency checks for T and
% dt
%
% @change{0,5,dw,2011-11-02} Modified the set.ODESolver method so that the MaxTimestep value is
% set to empty if implicit solvers are used. If again an explicit solver is used, a warning is
% issued if the models.BaseDynSystem.MaxTimestep value of the corresponding System is empty.
%
% @change{0,5,dw,2011-10-14} Removed the TimeDirty flag as it wasnt used properly/at all.
%
% @change{0,5,dw,2011-09-29}
% - New flag-field RealTimePlotting that calls the new plotSingle
% method in order to display the system as it is simulated.
% - Made ODESolver dependent to implement connection to
% RealTimePlotting.
% - New double field RealTimePlottingMinPause to enable timely display
% of in-simulation states.
%
% @change{0,3,sa,2011-05-10} Implemented setters for the rest of the
% properties
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
% 
% @change{0,3,dw,2011-04-15} Added a dependent GScaled property that returns the norm-inducing
% matrix G scaled with the current System.StateScaling property.
%
% @change{0,3,dw,2011-04-05} 
% - Removed the getConfigStr-Method and moved it to
% Utils.getObjectConfig
% - Added a setter for System checking for self-references
%
% @new{0,2,dw,2011-03-08} Implemented time scaling via addition of the
% property models.BaseModel.tau and dependent attributes
% models.BaseModel.dtscaled and models.BaseModel.Tscaled. This
% way model data can be entered in original units and the
% system calculates with the scaled time values. The main change is in
% @ref models.BaseModel.computeTrajectory where the ODE solver is
% called with the scaled time steps and the resulting timesteps are
% re-scaled to their original unit.
%
% @change{0,1,dw} Generalized scalar product via `<x,y>_G = x^tGy`,
% default `I_d` for `d\in\N`
%
% @new{0,1,dw} String output of all model settings via method
% getObjectConfig
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing    
    
    properties(SetObservable)
        % The actual dynamical system used in the model.
        %
        % @propclass{critical} No simulations without dynamical system.
        %
        % @default []
        System = [];
        
        % The name of the Model
        %
        % @propclass{optional}
        Name = 'Base Model';
                              
        % The custom scalar product matrix `\vG`
        %
        % In some settings the state variables have a special meaning (like
        % DOF's in FEM simulations) where the pure `L^2`-norm has less
        % meaning than a custom norm induced by a symmetric positive
        % definite matrix G. If `d\in\N` is the number of state
        % variables (i.e. dimensions of `\vx(t)`), then we must have
        % `\vG\in\R^{d\times d}`.
        %
        % Leave at default value `1` if `G=I_d` should be assumed.
        %
        % @propclass{optional}
        %
        % @default 1 @type matrix<double>
        G = 1;
        
        % Minimum pause between successive steps when RealTimePlotting is
        % enabled.
        %
        % @propclass{optional} Changes pause length between timesteps
        %
        % @default .1 @type double
        %
        % @see RealTimePlotting
        RealTimePlottingMinPause = .1;
    end
    
    properties(SetAccess=private, Dependent)
        % Evaluation points `\{0=t_0,\ldots,t_n=T\}` of the model 
        Times;
        
        % The time steps Times in scaled time units `\tilde{t_i} = \frac{t_i}{\tau}`
        %
        % See also: tau
        %
        % @default Times @type rowvec<double>
        scaledTimes;
        
        % The scaled end time `\tilde{T} = \frac{T}{\tau}`
        %
        % See also: tau T
        %
        % @default T @type double
        Tscaled;
    end
    
    properties(Dependent, SetObservable)
        % Time scaling factor `\tau`
        %
        % If used, the values from T and dt are getting scaled by tau when
        % calling simulate.
        %
        % @propclass{scaling}
        %
        % @default 1 @type double
        tau;
        
        % The final timestep `T` up to which to simulate.
        %
        % NOTE: When changing this property any offline computations have
        % to be repeated in order to obtain a new reduced model.
        %
        % @propclass{important} Defines the end time `T` up to which the dynamical system has to be
        % simulated.
        %
        % @type double @default 1
        T;
        
        % The desired time-stepsize `\Delta t` for simulations.
        %
        % @attention - This property is influencing the resulting output-times at which the
        % dynamical system is computed. If you need to set a maximum time-step size due to CFL
        % conditions, for example, use the models.BaseDynSystem.MaxTimeStep property.
        % - When changing this property any offline computations have
        % to be repeated in order to obtain a new reduced model.
        %
        % @propclass{critical}
        %
        % @default 0.1 @type double
        %
        % See also: dtscaled
        dt;
        
        % The solver to use for the ODE.
        % Must be an instance of any solvers.BaseSolver subclass.
        %
        % See also: solvers BaseSolver ode23 ode45 ode113
        %
        % @propclass{important} Choose an appropriate ODE solver for your
        % system.
        %
        % @type solvers.BaseSolver @default solvers.MLWrapper(ode23)
        ODESolver;
        
        % Determines if the simulation should plot intermediate steps
        % during computation.
        %
        % Disabled by default.
        %
        % @propclass{optional} Additionally displays the system's plot
        % during simulations.
        %
        % @default false @type logical
        RealTimePlotting;
    end
    
    properties(SetAccess=private)
        % The scaled timestep `\tilde{\Delta t} = \frac{\Delta t}{\tau}`
        %
        % @note Due to performance reasons this property is not computed
        % dependently but fitted any time dt or tau are changed.
        %
        % If tau is used, this value is `\tilde{dt} = dt/\tau`
        %
        % @default .1 (as in dt)
        % See also: tau dt
        dtscaled = .1;
        
        % Contains the GIT revision of this model when it was last saved.
        %
        % This field is automatically set upon saveobj.
        %
        % @type char @default ''
        gitRefOnSave = '';
    end
    
    properties(Access=private)
        ftau = 1;
        fdt = .1;
        fT = 1;
        frtp = false;
        fODEs;
        steplistener;
        ctime;
    end
    
    methods
        
        function this = BaseModel
            
            % Call superclass constructor first
            % (not necessary in this version as automatically called first, but one never knows..)
            this = this@KerMorObject;
            
            % Init defaults
            this.ODESolver = solvers.MLWrapper(@ode23);
            
            % Register default properties
            this.registerProps('System','T','ODESolver','dt','G',...
                'tau','RealTimePlottingMinPause','RealTimePlotting');
        end
        
        function [t, y, sec, x] = simulate(this, mu, inputidx)
            % Simulates the system and produces the system's output.
            %
            % Both parameters are optional. (Which to provide will be
            % determined by the actual system anyways)
            %
            % Parameters:
            % mu: The concrete mu parameter sample to simulate for.
            % inputidx: The index of the input function to use.
            %
            % Return values:
            % t: The times at which the model was evaluated
            % y: Depending on the existance of an output converter, this
            %    either returns the full trajectory or the processed output
            %    at times t.
            % sec: the seconds needed for simulation.
            % x: The (scaled, if used) state space variables before output conversion.
            %
            % @todo: fix Input checks (set inidx=1 iff one input is
            % there, otherwise error)
            % @todo: switch return arguments sec & x + tests
            if nargin < 3
                inputidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            this.WorkspaceVariableName = inputname(1);
            
            if this.RealTimePlotting
                this.ctime = tic;
            end
            % Get scaled state trajectory
            [t, x, time] = this.computeTrajectory(mu, inputidx);
            
            % Measure rest of time
            starttime = tic;
            
            % Convert to output
            y = this.System.computeOutput(x);
            
            % Scale times back to real units
            t = t*this.tau;
            
            sec = toc(starttime) + time;
            
            this.WorkspaceVariableName = '';
        end
        
        function [f, ax] = plot(this, t, y, varargin)
            % Plots the results of the simulation.
            % Override in subclasses for a more specific plot if desired.
            %
            % Parameters:
            % t: The simulation times `t_i` @type rowvec
            % y: The simulation output matrix `y`, i.e. `y(t_i)` @type matrix
            % varargin: Any further arguments for customized plots
            %
            % Return values:
            % f: The figure handle @type handle
            % ax: The axes handle @type handle
            if isempty(varargin)
                f = figure;
                ax = gca(f);
            else
                ax = varargin{1};
                f = get(ax,'Parent');
            end
            y = Utils.preparePlainPlot(y);
            plot(ax,t,y);
            title(ax,sprintf('Output plot for model "%s"', this.Name));
            xlabel(ax,'Time'); ylabel(ax,'Output values');
        end
        
        function [f, ax] = plotState(this, t, x, varargin)
            % Plots the results of the simulation.
            % Override in subclasses for a more specific plot if desired.
            %
            % Parameters:
            % t: The simulation times `t_i` @type rowvec
            % y: The simulation output matrix `y`, i.e. `y(t_i)` @type matrix
            % varargin: Any further arguments for customized plots
            %
            % Return values:
            % f: The figure handle @type handle
            % ax: The axes handle @type handle
            if isempty(varargin)
                f = figure;
                ax = gca(f);
            else
                ax = varargin{1};
                f = get(ax,'Parent');
            end
            x = Utils.preparePlainPlot(x);
            plot(ax,t,x);
            title(ax,sprintf('State space plot for model "%s"', this.Name));
            xlabel(ax,'Time'); ylabel(ax,'State space values');
        end
        
        function plotSingle(this, t, y, varargin)
            % Plots a single solution.
            % Override in subclasses for specific plot behaviour.
            %
            % The default method is simply to use the full plot default
            % method.
            %
            % Parameters:
            % t: The current time `t`
            % y: The system's output `y(t)`
            % varargin: Any arguments that should be passed on to inner plotting methods (model
            % dependent)
            this.plot(t, y, varargin{:});
        end
        
        function [t, x, ctime] = computeTrajectory(this, mu, inputidx)
            % Computes a solution/trajectory for the given mu and inputidx in the SCALED state
            % space.
            %
            % Parameters:
            % mu: The parameter `\mu` for the simulation
            % inputidx: The integer index of the input function to use. If
            % more than one inputs are specified this is a necessary
            % argument.
            %
            % Return values:
            % t: The times at which the model was evaluated. Will equal the property Times
            % @type rowvec
            % x: The state variables at the corresponding times t. @type matrix<double>
            % ctime: The time needed for computation. @type double
            %
            % @change{0,7,dw,2013-03-19} Fixed the MaxStep setting of the ODE solver when
            % time-scaling is used.
            if nargin < 3
                inputidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            sys = this.System;
            
            % Stop the time
            st = tic;
            
            % Prepare the system by setting mu and inputindex.
            sys.setConfig(mu, inputidx);
            
            %% Solve ODE
            slv = this.ODESolver;
            slv.MaxStep = []; slv.InitialStep = [];
            if ~isempty(sys.MaxTimestep)
                % Remember: When scaling is used, these are the 
                slv.MaxStep = sys.MaxTimestep/this.tau;
                slv.InitialStep = .5*sys.MaxTimestep/this.tau;
            end
            
            % Assign jacobian information if available and solver is
            % implicit
            if isa(slv,'solvers.AImplSolver')
                % Set jacobian
                slv.JPattern = [];
                slv.JacFun = [];
                if ~isempty(sys.A) && ~isempty(sys.f)
                    slv.JacFun = @(t, x)sys.A.getStateJacobian(x, t, mu) + sys.f.getStateJacobian(x, t, mu);
                    if ~isempty(sys.A.JSparsityPattern) && ~isempty(sys.f.JSparsityPattern)
                        [i,j] = find(sys.A.JSparsityPattern + sys.f.JSparsityPattern);
                            slv.JPattern = sparse(i,j,ones(length(i),1),...
                                sys.A.fDim,sys.A.xDim);
                    end
                elseif ~isempty(sys.A)
                    slv.JacFun = @(t, x)sys.A.getStateJacobian(x, t, mu);
                    if  ~isempty(sys.A.JSparsityPattern)
                        slv.JPattern = sys.A.JSparsityPattern;
                    end
                elseif ~isempty(sys.f)
                    slv.JacFun = @(t, x)sys.f.getStateJacobian(x, t, mu);
                    if  ~isempty(sys.f.JSparsityPattern)
                        slv.JPattern = sys.f.JSparsityPattern;
                    end
                end
            end
            
            % Assign mass matrix to solver if present
            slv.M = [];
            if ~isempty(sys.M)
                slv.M = sys.M;
            end
            
            % Determine correct ODE function (A,f,B combination)
            str = '';
            if ~isempty(sys.A)
                str = [str 'A'];
            end
            if ~isempty(sys.f)
                str = [str 'f'];
            end
            if ~isempty(sys.B) &&  ~isempty(inputidx)
                str = [str 'B'];
            end
            odefun = eval(sprintf('@(t,x)sys.ODEFun_%s(t,x)',str));
            
            % Call solver
            [t, x] = slv.solve(odefun, this.scaledTimes, this.getX0(mu));
            
            % Get used time
            ctime = toc(st);
        end
        
        function mu = getRandomParam(this, num, seed)
            % Gets a random parameter sample from the system's parameter
            % domain P
            %
            % Parameters:
            % num: The number of random parameters to return. @type integer @default 1
            % seed: The seed for the number generator. Leave empty for cputime initialization.
            % @type integer @default 'round(cputime*100)'
            %
            % Return values:
            % mu: A matrix of random parameters within the specified range for this model.
            % @type matrix<double>
            %
            % @change{0,6,dw,2012-07-13} Moved this method from models.BaseDynSystem to here
            % for more convenience.
            if nargin < 3
                seed = round(cputime*100);
                if nargin < 2
                    num = 1;
                end
            end
            s = this.System;
            r = RandStream('mt19937ar','Seed',seed);
            if s.ParamCount > 0
                pmin = [s.Params(:).MinVal]';
                pmax = [s.Params(:).MaxVal]';
                mu = r.rand(s.ParamCount,num) .* repmat(pmax-pmin,1,num) + repmat(pmin,1,num);
            else
                mu = [];
            end
        end
    end
    
    methods(Access=protected)
        function x0 = getX0(this, mu)
            % Gets the initial state variable at `t=0`.
            %
            % This is exported into an extra function as it gets overridden
            % in the ReducedModel subclass, where ErrorEstimators possibly
            % change the x0 dimension.
            %
            % Parameters:
            % mu: The parameter `\mu` to evaluate `x_0(\mu)`. Use [] for
            % none.
            ss = this.System.StateScaling;
            if size(mu,2) > 1
                ss = repmat(this.System.StateScaling,1,size(mu,2));
            end
            x0 = this.System.x0.evaluate(mu) ./ ss;
        end
        
        function this = saveobj(this)
            % Store the current GIT branch in the object.
            this.gitRefOnSave = KerMor.getGitBranch;
        end
    end
    
    methods(Access=protected, Sealed)
        function plotstep(this, src, ed)%#ok
            % Callback for the ODE solvers StepPerformed event that enables
            % during-simulation-plotting.
            y = this.System.computeOutput(ed.States);
            this.plotSingle(ed.Times * this.tau,y);
            drawnow;
            % Gets first set in "simulate"
            per = toc(this.ctime);
            pause(max(0,this.RealTimePlottingMinPause-per));
            this.ctime = tic;
        end
    end
    
    %% Getter & Setter
    methods
        function value = get.Times(this)
            value = 0:this.dt:this.T;
        end
        
        function value = get.scaledTimes(this)
            value = 0:this.dtscaled:this.Tscaled;
        end
        
        function value = get.Tscaled(this)
            if isempty(this.tau)
                value = this.T;
            else
                value = this.T/this.tau;
            end
        end
        
        function dt = get.dt(this)
            dt = this.fdt;
        end
        
        function tau = get.tau(this)
            tau = this.ftau;
        end
        
        function set.T(this, value)
            if ~isscalar(value) || value < 0
                error('T must be a positive real scalar.');
            elseif value < this.fdt
                error('Timestep dt must be smaller or equal to T.');
            end
            if value ~= this.fT
                this.fT = value;
            end
        end
        
        function set.dt(this, value)
            if ~isscalar(value) || value <= 0
                error('dt must be a positive real scalar.');
            elseif value > this.fT
                error('Timestep dt must be smaller or equal to T.');
            end
            if this.fdt ~= value
                this.dtscaled = value/this.ftau;
                this.fdt = value;
            end
        end
        
        function set.tau(this, value)
            if ~isreal(value) ||~isscalar(value)
                error('tau must be a positive real scalar.');
            end
            if this.ftau ~= value
                this.dtscaled = this.fdt/value;
                this.ftau = value;
            end
        end
        
        function set.System(this, value)
            % @note Usually an empty system is not allowed. But as this is a superclass for
            % both full and reduced models, one
            if ~isempty(value) && ~isa(value,'models.BaseDynSystem')
                error('The System property must be a subclass of models.BaseDynSystem.');
            end
            if (isequal(this,value))
                warning('KerMor:selfReference','Careful with self-referencing classes. See BaseFullModel class documentation for details.');
            end
            this.System = value;
        end
        
        function set.G(this, value)
            % @todo check for p.d. and symmetric, -> sparsity?
            if isempty(value)
                error('G must not be empty. Use G=1 for default euclidean scalar product and norms');
            elseif any(abs(value-value') > eps)
                error('G must be symmetric.');
            end
            this.G = value;
        end
        
        function set.ODESolver(this, value)
            this.checkType(value,'solvers.BaseSolver');
            % Add listener if new ODE solver is passed and real time
            % plotting is turned on.
            if this.frtp && this.fODEs ~= value
                if ~isempty(this.steplistener)
                    delete(this.steplistener);
                end
                this.steplistener = value.addlistener('StepPerformed',@this.plotstep);
            end
            % Disable the MaxTimestep value if an implicit solver is used.
            if ~isempty(this.System) 
                if isa(value,'solvers.AImplSolver')
                    if ~isempty(this.System.MaxTimestep)
                        fprintf('BaseModel: Disabling system''s MaxTimestep due to use of an implicit solver.\n');
                    end
                    this.System.MaxTimestep = [];
                elseif isempty(this.System.MaxTimestep)
                    warning('KerMor:BaseModel','Attention: Setting an explicit solver for a system without MaxTimestep set. Please check.');
                end
            end
            this.fODEs = value;
        end
        
        function value = get.ODESolver(this)
            value = this.fODEs;
        end
        
        function set.Name(this, value)
            if ~ischar(value)
                error('name is acharacter field');
            end
            this.Name = value;
        end
        
        function value = get.RealTimePlotting(this)
            value = this.frtp;
        end
        
        function value = get.T(this)
            value = this.fT;
        end
        
        function set.RealTimePlotting(this, value)
            if ~islogical(value)
                error('Value must be boolean.');
            end
            if value
                if isempty(this.steplistener)
                    this.steplistener = this.ODESolver.addlistener('StepPerformed',@this.plotstep);
                end
            else
                if ~isempty(this.steplistener)
                    delete(this.steplistener);
                    this.steplistener = [];
                end
            end
            this.frtp = value;
        end
    end
end

