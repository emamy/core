classdef BaseModel < handle
    % Base class for both full and reduced models.
    %
    % This class gathers all common functionalities of models in the
    % KerMor framework.
    % The most important method would be @code [t,y] =
    % simulate(mu,inputidx) @endcode which computes the system's solution
    % for given `\mu` and input number (if applicable).  Also a plot
    % wrapper is provided that refers to the plotting methods within the
    % model's system.
    % Test: `x\in\R`
    %
    % @author Daniel Wirtz @date 19.03.2010
    %
    % @new{0,3,2011-03-08} Implemented time scaling via addition of the
    % property @ref models.BaseModel.tau "tau" and dependent attributes @ref models.BaseModel.dtscaled "dtscaled" and @ref
    % models.BaseModel.Tscaled "Tscaled". This way model data can be entered in original units and the
    % system calculates with the scaled time values. The main change is in
    % @ref models.BaseModel.computeTrajectory where the ODE solver is
    % called with the scaled time steps and the resulting timesteps are
    % re-scaled to their original unit.
    %
    % @change{0,1} Generalized scalar product via `<x,y>_G = x^tGy`,
    % default `I_d` for `d\in\N`
    % @change{0,1} String output of all model settings via method
    % getObjectConfig
    
    properties
        % The actual dynamical system used in the model.
        System;
        
        % The name of the Model
        Name = 'Base Model';
        
        % The verbose output level at simulations
        Verbose = 0;
        
        % The final timestep up to which to simulate.
        %
        % NOTE: When changing this property any offline computations have
        % to be repeated in order to obtain a new reduced model.
        T = 1;
        
        % The desired time-stepsize for simulations
        %
        % NOTE: When changing this property any offline computations have
        % to be repeated in order to obtain a new reduced model.
        %
        % @default 0.1
        % See also: dtscaled
        dt = .1;
        
        % The solver to use for the ODE.
        % Must be an instance of any solvers.BaseSolver subclass.
        %
        % Default: @code solvers.MLWrapper(@ode23) @endcode
        %
        % See also: solvers BaseSolver ode23 ode45 ode113
        ODESolver;
        
        % The custom scalar product matrix
        %
        % In some settings the state variables have a special meaning (like
        % DOF's in FEM simulations) where the pure `L^2`-norm has less
        % meaning than a custom norm induced by a symmetric positive
        % definite matrix G. If `d\in\mathbb{N}` is the number of state
        % variables (i.e. dimensions of `x(t)`), then we must have
        % `G\in\mathbb{R}^{d\times d}`.
        %
        % Leave at `1` if `G=I_d` should be assumed.
        %
        % Default:
        % 1
        G = 1;
        
        % Time scaling factor
        %
        % If used, the values from T and dt are getting scaled by tau when
        % calling simulate.
        %
        % @default 1
        tau = 1;
    end
    
    properties(Dependent)
        % Evaluation points of the model as increasing array
        Times;
        
        % The time steps in scaled time units
        %
        % See also: tau
        scaledTimes;
        
        % The scaled end time T
        %
        % If tau is used, this value is `\tilde{T} = T/\tau`
        % See also: tau T
        Tscaled;
    end
    
    properties(Access=protected)
        % Flag that indicates changes in either T or dt after
        % offlineGenerations have been performed.
        TimeDirty;
    end
    
    properties(SetAccess=private)
        % The scaled timestep dt
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
    
    methods
        
        function this = BaseModel
             this.System = models.BaseDynSystem;
             this.ODESolver = solvers.MLWrapper(@ode23);
        end
        
        function [t,y,sec,x] = simulate(this, mu, inputidx)
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
                if ~isempty(this.System.Inputs)
                    warning('BaseModel:NoInputSpecified',['You must specify '...
                        'an input index if inputs are set up. Using inputidx=1']);
                    inputidx = 1;
                else
                    inputidx = [];
                end
                if nargin < 2
                    mu = [];
                end
            end
            if isempty(mu) && ~isempty(this.System.Params)
                error('A model with parameters cannot be simulated without a parameter argument.');
            end
            
            starttime = tic;
            
            % Get state trajectory
            [t,x] = this.computeTrajectory(mu, inputidx);
            
            % Convert to output
            y = this.System.C.computeOutput(t,x,mu);
            
            sec = toc(starttime);
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
        
        function cfg = getConfigStr(this)
            % Returns a string that contains all available class and
            % subclass information about the current model.
            
            %res = ['Configuration of model ''' this.Name '''' this.getObjectConfig(this, 0, 20)];
            %cfg = res;
            %cfg = sprintf(res);
            cfg = sprintf(['Configuration of model ''' this.Name '''' this.getObjectConfig(this, 0, 20)]);
        end
        
        function [t,x] = computeTrajectory(this, mu, inputidx)
            % Computes a solution/trajectory for the given mu and inputidx.
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
            % x: Depending on the existance of an output converter, this
            %    either returns the full trajectory or the processed output
            %    at times t.
            
            % Validity checks
            if this.System.InputCount > 0 && isempty(inputidx)
                if this.System.InputCount == 1
                    inputidx = 1;
                else
                    error('If more than one input u is set you must specify an inputidx.');
                end
            end
            
            % Setup simulation-time constant data (if available)
            if isa(this.System,'ISimConstants')
                this.System.updateSimConstants;
            end
            if isa(this.System.f,'ISimConstants')
                this.System.f.updateSimConstants;
            end

            % Get target ODE function
            if isempty(this.System.f)
                error('No system''s core function specified. ODE function creation impossible; first set "f" property.');
            end
            odefun = this.System.getODEFun(mu, inputidx);
            
            % Get initial x0
            x0 = this.getX0(mu);
            
            % Solve ODE
            if ~isempty(this.System.MaxTimestep)
                % Remember: When scaling is used, these are the 
                this.ODESolver.MaxStep = this.System.MaxTimestep;
                this.ODESolver.InitialStep = .5*this.System.MaxTimestep;
            end
            [t,x] = this.ODESolver.solve(odefun, this.scaledTimes, x0);
            % Scale times back to real units
            t = t*this.tau;
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
            x0 = this.System.x0(mu);
        end
    end
    
    methods(Access=protected,Sealed)
               
        function checkType(this, obj, type)%#ok
            % Object typechecker.
            % Checks if a given object is of the specified type and throws
            % an error if not.
            % Convenience method.
            if ~isempty(obj) && ~isa(obj, type)
                error(['Wrong type ''' class(obj) ''' for this property. Has to be a ' type]);
            end
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
        
        function set.T(this, value)
            if ~isposrealscalar(value)
                error('T must be a positive real scalar.');
            end
            if this.T ~= value
                this.T = value;
                this.TimeDirty = true;%#ok
            end
        end
        
        function set.dt(this, value)
            if ~isposrealscalar(value)
                error('T must be a positive real scalar.');
            end
            if this.dt ~= value
                if isempty(this.tau)%#ok
                    this.dtscaled = value;%#ok
                else
                    this.dtscaled = value/this.tau;%#ok
                end
                this.dt = value;
                this.TimeDirty = true;%#ok
            end
            
        end
        
        function set.tau(this, value)
            if ~isposrealscalar(value)
                error('tau must be a positive real scalar.');
            end
            if this.tau ~= value
                if ~isempty(this.dt)%#ok
                    this.dtscaled = this.dt/value;%#ok
                end
                this.tau = value;
                this.TimeDirty = true;%#ok
            end
        end
        
        function set.G(this, value)
            % @todo check for p.d. and symmetric, -> sparsity?
            this.G = value;
        end
        
        function set.ODESolver(this, value)
            this.checkType(value,'solvers.BaseSolver');%#ok
            this.ODESolver = value;
        end
    end
    
    methods(Access=private)
        function str = getObjectConfig(this, obj, numtabs, depth)
            str = '';
            if depth == 0
                warning('Un:important','Exceeded the maximal recursion depth. Ignoring deeper objects.');
                return;
            end
            mc = metaclass(obj);
            if isfield(obj,'Name')
                name = obj.Name;
            else
                name = mc.Name;
            end
            str = [str ': ' name '\n'];
            for idx = 1:length(mc.Properties)
                p = mc.Properties{idx};
                if strcmp(p.GetAccess,'public')
                    str = [str repmat('\t',1,numtabs)];%#ok
                    pobj = obj.(p.Name);
                    if ~isempty(pobj) && ~isequal(obj,pobj)
                        if isobject(pobj)
                            str = [str p.Name this.getObjectConfig(pobj, numtabs+1, depth-1)];%#ok
                        elseif isnumeric(pobj)
                            if numel(pobj) > 20
                                str = [str p.Name ': [' num2str(size(pobj)) '] ' class(pobj)];%#ok
                            else
                                pobj = reshape(pobj,1,[]);
                                str = [str p.Name ': ' num2str(pobj)];%#ok
                            end
                        elseif isstruct(pobj)
                            if any(size(pobj) > 1)
                                str = [str p.Name ' (struct, [' num2str(size(pobj)) ']), fields: '];%#ok
                            else
                                str = [str p.Name ' (struct), fields: '];%#ok
                            end
                            fn = fieldnames(pobj);
                            for fnidx = 1:length(fn)
                                %str = [str fn{fnidx} this.getObjectConfig(pobj, numtabs+1, depth-1)];%#ok
                                str = [str fn{fnidx} ','];%#ok
                            end
                            str = str(1:end-1);
                        elseif isa(pobj,'function_handle')
                            str = [str p.Name ' (f-handle): ' func2str(pobj)];%#ok
                        end
                    else
                        if isequal(obj,pobj)
                            str = [str p.Name ': self-reference'];%#ok
                        else
                            str = [str p.Name ': empty'];%#ok
                        end
                    end
                    str = [str '\n'];%#ok
                    %str = strcat(str,'\n');
                end
            end
            % Take away the last \n
            str = str(1:end-2);
        end
    end
end
