classdef BaseModel < KerMorObject
    % Base class for both full and reduced models.
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
        
        % The solver to use for the ODE.
        % Must be an instance of any solvers.ode.BaseSolver subclass.
        %
        % Default: @code solvers.ode.MLWrapper(@ode23) @endcode
        %
        % See also: solvers BaseSolver ode23 ode45 ode113
        %
        % @propclass{important}
        ODESolver;
        
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
    end
    
    properties(Access=protected)
        % Flag that indicates changes in either T or dt after
        % offlineGenerations have been performed.
        %
        % @todo check if this property still makes sense
        TimeDirty;
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
    end
    
    methods
        
        function this = BaseModel
            
            % Call superclass constructor first
            % (not necessary in this version as automatically called first, but one never knows..)
            this = this@KerMorObject;
            
            % Init defaults
            this.ODESolver = solvers.ode.MLWrapper(@ode23);
            
            % Register default properties
            this.registerProps('System','T','ODESolver','dt','G','tau');
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
            
            % Get scaled state trajectory
            [t,x] = this.computeTrajectory(mu, inputidx);
            
            % Scale times back to real units
            x = bsxfun(@times,x,this.System.StateScaling);
            t = t*this.tau;
            
            % Convert to output
            y = this.System.C.computeOutput(t,x,mu);
            
            sec = toc(starttime);
            
            this.WorkspaceVariableName = '';
        end
        
        function plot(this, t, y, ax)
            % Plots the results of the simulation.
            % Override in subclasses for a more specific plot if desired.
            %
            % Parameters:
            % t: The simulation times `t_i`
            % y: The simulation output matrix `y`, i.e. `y(t_i)`
            % ax: An axes handle as plotting target. Optional; if not
            % specified a new figure is created.
            if nargin < 4
                figure;
                ax = gca;
            end
            y = general.Utils.preparePlainPlot(y);
            plot(ax,t,y);
            title(ax,sprintf('Plot for output of model "%s"', this.Name));
            xlabel(ax,'Time'); ylabel(ax,'Output functions');
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
            
            %% Setup simulation-time constant data (if available)
            if isa(this.System,'ISimConstants')
                this.System.prepareConstants;
            end
            if isa(this.System.f,'ISimConstants')
                this.System.f.prepareConstants;
            end
  
            %% Pass mu and input to system
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
%             target.TimeDirty = this.TimeDirty;
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
            if ~isposrealscalar(value)
                error('T must be a positive real scalar.');
            end
            if this.T ~= value
                this.T = value;
                this.TimeDirty = true;%#ok -> setter does not change anything
            end
        end
        
        function set.dt(this, value)
            if ~isposrealscalar(value)
                error('dt must be a positive real scalar.');
            end
            if this.fdt ~= value
                this.dtscaled = value/this.ftau;
                this.fdt = value;
                this.TimeDirty = true;
            end
        end
        
        function set.tau(this, value)
            if ~isposrealscalar(value)
                error('tau must be a positive real scalar.');
            end
            if this.ftau ~= value
                this.dtscaled = this.fdt/value;
                this.ftau = value;
                this.TimeDirty = true;
            end
        end
        
        function set.System(this, value)
            if ~isa(value,'models.BaseDynSystem')
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
            this.checkType(value,'solvers.ode.BaseSolver');%#ok
            this.ODESolver = value;
        end
        
        function set.Name(this, value)
            if ~ischar(value)
                error('name is acharacter field');
            end
            this.Name = value;
        end
        
        function set.TimeDirty(this, value)
            if ~islogical(value)
                error('value must be a logical');
            end
            this.TimeDirty = value;
        end
        
    end
end

