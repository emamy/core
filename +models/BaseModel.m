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
% general.Utils.getObjectConfig
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
                
        % The final timestep `T` up to which to simulate.
        %
        % NOTE: When changing this property any offline computations have
        % to be repeated in order to obtain a new reduced model.
        %
        % @propclass{important} Defines the end time `T` up to which the dynamical system has to be
        % simulated.
        T = 1;
              
        % The custom scalar product matrix `G`
        %
        % In some settings the state variables have a special meaning (like
        % DOF's in FEM simulations) where the pure `L^2`-norm has less
        % meaning than a custom norm induced by a symmetric positive
        % definite matrix G. If `d\in\mathbb{N}` is the number of state
        % variables (i.e. dimensions of `x(t)`), then we must have
        % `G\in\mathbb{R}^{d\times d}`.
        %
        % Leave at default value `1` if `G=I_d` should be assumed.
        %
        % @propclass{optional}
        %
        % @default 1
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
        scaledTimes;
        
        % The scaled end time `\tilde{T} = \frac{T}{\tau}`
        %
        % See also: tau T
        Tscaled;
        
        % The scaled version of G.
        %
        % Equals `\tilde{G} = D'GD` for `D=diag(s)` and `s` being the System.StateScaling property.
        %
        % Use this whenever having to take the real G-norm of some scaled state variables
        %
        % See also: G System.StateScaling
        GScaled;
    end
    
    properties(Dependent, SetObservable)
        % Time scaling factor `\tau`
        %
        % If used, the values from T and dt are getting scaled by tau when
        % calling simulate.
        %
        % @propclass{scaling}
        %
        % @default 1
        tau;
        
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
        % @default 0.1
        %
        % See also: dtscaled
        dt;
        
        % The solver to use for the ODE.
        % Must be an instance of any solvers.ode.BaseSolver subclass.
        %
        % Default: @code solvers.ode.MLWrapper(@ode23) @endcode
        %
        % See also: solvers BaseSolver ode23 ode45 ode113
        %
        % @propclass{important} Choose an appropriate ODE solver for your
        % system.
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
    
    properties(SetAccess=private, GetAccess=protected)
        % The variable name of the variable denoting the current model instance.
        %
        % Will be set upon calling 'simulate' and will equal '' outside the scope of a simulation.
        %
        % Can be used in subclasses to create useful links regarding functions of the model, i.e. is
        % used in BaseFullModels checkProperties member to create a link for the critical properties
        % summary.
        %
        % @default ''
        WorkspaceVariableName = '';
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
    end
    
    properties(Access=private)
        ftau = 1;
        fdt = .1;
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
            this.ODESolver = solvers.ode.MLWrapper(@ode23);
            
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
            % x: The state space variables before output conversion.
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
            
            starttime = tic;
            
            if this.RealTimePlotting
                this.ctime = tic;
            end
            % Get scaled state trajectory
            [t,x] = this.computeTrajectory(mu, inputidx);
            % Re-scale state variable
            x = bsxfun(@times,x,this.System.StateScaling);
            % Convert to output
            y = this.System.C.computeOutput(t,x,mu);
            
            % Scale times back to real units
            t = t*this.tau;
            
            sec = toc(starttime);
            
            this.WorkspaceVariableName = '';
        end
        
        function plot(this, t, y)
            % Plots the results of the simulation.
            % Override in subclasses for a more specific plot if desired.
            %
            % Parameters:
            % t: The simulation times `t_i`
            % y: The simulation output matrix `y`, i.e. `y(t_i)`
            y = general.Utils.preparePlainPlot(y);
            plot(ax,t,y);
            title(ax,sprintf('Plot for output of model "%s"', this.Name));
            xlabel(ax,'Time'); ylabel(ax,'Output functions');
        end
        
        function plotSingle(this, t, y, ax)
            % Plots a single solution.
            % Override in subclasses for specific plot behaviour.
            %
            % The default method is simply to use the full plot default
            % method.
            %
            % Parameters:
            % t: The current time `t`
            % y: The system's output `y(t)`
            % ax: The target axes to plot into. @default []
            if nargin < 4
                this.plot(t, y);
            else
                this.plot(t, y, ax);
            end
        end
        
        function [t, x] = computeTrajectory(this, mu, inputidx)
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
            % t: The times at which the model was evaluated. Will equal
            % the property Times
            % x: The state variables at the corresponding times t.
            if nargin < 3
                inputidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            
            % Prepare the system by setting mu and inputindex.
            this.System.setConfig(mu, inputidx);
            
            %% Solve ODE
            slv = this.ODESolver;
            if ~isempty(this.System.MaxTimestep)
                % Remember: When scaling is used, these are the 
                slv.MaxStep = this.System.MaxTimestep;
                slv.InitialStep = .5*this.System.MaxTimestep;
            end
            % Assign jacobian evaluation function if available
            if isa(slv,'solvers.ode.BaseImplSolver')
                if isa(this.System.f,'dscomponents.IJacobian')
                    slv.JacFun = @(t, x)this.System.f.getStateJacobian(x, t, mu);
                else
                    slv.JacFun = [];
                end
            end
            % Call solver
            [t, x] = slv.solve(@(t,x)this.System.ODEFun(t,x), this.scaledTimes, this.getX0(mu));
        end
        
%         function target = clone(this, target)            
%             % Clones this instance into another instance given by target.
%             if nargin < 2 || ~isa(target, 'models.BaseModel')
%                 error('The target argument must be given and a valid BaseModel subclass.');
%             end
%             % Clone the system & odesolver, too
%             target.System = this.System.clone;
%             target.ODESolver = this.ODESolver.clone;
%             target.Name = this.Name;
%             target.Verbose = this.Verbose;
%             target.T = this.T;
%             target.dt = this.dt;
%             target.G = this.G;
%             target.tau = this.tau;
%             target.dtscaled = this.dtscaled;
%         end
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
            x0 = this.System.x0.evaluate(mu) ./ this.System.StateScaling;
        end
    end
    
    methods(Access=protected, Sealed)
        function plotstep(this, src, ed)%#ok
            % Callback for the ODE solvers StepPerformed event that enables
            % during-simulation-plotting.
            y = this.System.C.computeOutput(ed.Times, ...
                this.System.StateScaling .* ed.States, this.System.mu);
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
        
        function gs = get.GScaled(this)
            gs = diag(this.System.StateScaling) * this.G * diag(this.System.StateScaling);
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
            end
            this.T = value;
        end
        
        function set.dt(this, value)
            if ~isscalar(value) || value <= 0
                error('dt must be a positive real scalar.');
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
            this.G = value;
        end
        
        function set.ODESolver(this, value)
            this.checkType(value,'solvers.ode.BaseSolver');
            % Add listener if new ODE solver is passed and real time
            % plotting is turned on.
            if this.frtp && this.fODES ~= value
                if ~isempty(this.steplistener)
                    delete(this.steplistener);
                end
                this.steplistener = value.addlistener('StepPerformed',@this.plotstep);
            end
            % Disable the MaxTimestep value if an implicit solver is used.
            if ~isempty(this.System) 
                if isa(value,'solvers.ode.BaseImplSolver')
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

