classdef BaseFirstOrderSystem < KerMorObject
    % Base class for all KerMor first-order dynamical systems.
    %
    % To setup custom dynamical systems, inherit from this class.
    %
    % No system components are set by default, except a scalar
    % dscomponents.LinearOutputConv instance for simple state-space to
    % output forwarding.
    %
    % Programming convention:
    % Any system-specific properties should be added in the subclass, even
    % though sometimes one might tend to put settings into the core
    % function f (e.g. for space-discretized systems and thus all
    % space-discretization relevant settings). In order to have all
    % settings combined at one class one should rather pass the custom
    % system's instance to the ACoreFun implementing classes' constructor
    % and read settings from that instead.
    % For an example see the models.pcd.PCDSystem2D and pcd.CoreFun2D.
    %
    % @author Daniel Wirtz @date 2010-03-17
    %
    % @new{0,6,dw,2011-12-06} Added a new component M, an optional system's
    % mass matrix `M(t,\mu)x'(t) = f\ldots`.
    %
    % @change{0,5,dw,2011-10-17}
    % - Removed the ISimConstants interface as the
    % concept is too confusing when reducing models (core function
    % evaluations must be fully self-contained regarding current mu and
    % especially inputs)
    % - Renamed the prepareConstants to BaseDynSystem.setConfig
    % - New BaseDynSystem.getParamInfo method to display formatted
    % information about a given parameter vector.
    %
    % @change{0,5,dw,2011-07-07} New method getParamIndexFromName and fixed setConfig checks.
    %
    % @change{0,4,dw,2011-05-29} Re-changed the `x0(\mu)` property back to be a function handle. Higher
    % flexibility during subsequent simulations are to be rated higher than proper setting of the
    % property, which should be assumed is being done (especially) for `x0(\mu)`
    %
    % @change{0,4,dw,2011-05-13} Created new transient properties BaseDynSystem.mu and
    % BaseDynSystem.u in order to store the current `\mu` and `u(t)` used for simulations on
    % a higher level. Removed the old 'getODEFun' function and replaced it by a direct
    % implementation 'ODEFun'. This enables to avoid nested function handles with in turn allow for
    % a speedup of reduced simulations by almost a factor of 2.
    %
    % @change{0,3,dw,2011-05-03} Moved the x0 property into an abstract function. This way
    % subclasses must implemement the x0 function and is not only warned by the DPCS that it should
    % be set.
    %
    % @change{0,3,dw,2011-04-6} This class inherits from KerMorObject now
    % in order to enable correct saving/loading
    %
    % @new{0,3,dw,2011-04-05} Added a setter for property BaseDynSystem.f
    % that checks for self references.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetObservable)
        % The core f function from the dynamical system.
        %
        % @propclass{important}
        %
        % @type dscomponents.ACoreFun @default []
        f;
        
        % Represents a linear or affine-linear component of the dynamical
        % system.
        %
        % If there is any linear part in your system, assign it here to
        % take advantage of more involved solvers like semi-implicit euler.
        %
        % @propclass{optional}
        %
        % @type dscomponents.AffLinCoreFun,dscomponents.LinearCoreFun
        % @default []
        A = [];
        
        % The input conversion
        %
        % @propclass{optional}
        %
        % @default [] @type dscomponents.AInputConv
        B = [];
        
        % The output conversion
        % Defaults to an LinearOutputConv instance using a 1-matrix, which
        % just forwards the state variables and supports projection.
        %
        % @propclass{optional}
        %
        % See also: dscomponents.LinearOutputConv.project
        %
        % @type dscomponents.LinearOutputConv @default 1
        C;
        
        % Function handle to initial state evaluation.
        % 
        % The handle returns the initial value `x(0) = x_0` for the ODE
        % solver to start simulations from.
        %
        % Function parameter is the current parameter `\mu` and return
        % value is initial system state `x_0`.
        %
        % @propclass{critical} The initial value greatly influences the
        % simulation results.
        %
        % @type dscomponents.AInitialValue
        % @default []
        x0 = [];
        
        % The system's mass matrix.
        %
        % @default [] @type dscomponents.AMassMatrix
        M = [];
        
        % The system's algebraic constraints function
        %
        % @default [] @type dscomponents.ACoreFun
        g = [];
        
        % The system's possible input functions.
        % A cell array of function handles, each taking a time argument t.
        %
        % Note:
        % Adding an input to a system does not mean that any associated full model uses them in the
        % offline phase. In order to use an input during offline phase add the corresponding index
        % to the models.BaseFullModel.TrainingInputs property.
        %
        % @propclass{optional}
        Inputs = {};
        
        % The parameters usable for the dynamical system.
        %
        % @propclass{optional}
        Params = data.ModelParam.empty;
        
        % Indices of the parameter vector that are effectively used in the system's core
        % function.
        %
        % Set this property to the list of indices whose values are actually used at the
        % evaluation of the `f` function. Leaving it empty means none are used, and -1 means
        % all values are used (default).
        %
        % When building kernel expansion approximations this setting will be used to configure
        % the parameter kernel of the kernels.ParamTimeKernelExpansion.
        %
        % @propclass{critical} Using more `\mu` parameter vector elements in the kernel
        % approximations than actually required in the core function `f` introduces an
        % additional dependency of the nonlinearity on extra parameters which is not given in
        % the full model's core function.
        %
        % @type integer @default -1
        %
        % See also: kernels.ParamTimeKernelExpansion
        DependentParamIndices = -1;
        
        % The maximum timestep allowed for any ODE solvers.
        %
        % This might be necessary if the Core function encapsulates a
        % spatial discretization and thus CFL conditions apply, for example.
        %
        % @attention ''IMPORTANT!'' If any time scaling is used within the
        % model, this value must correspond to the ''scaled'' maximum
        % timestep.
        %
        % @propclass{critical} Too large time-steps might violate e.g. CFL
        % conditions and render the solution false and useless.
        %
        % @type double
        %
        % See also: BaseModel BaseModel.tau BaseModel.dtscaled
        MaxTimestep = [];
        
        % The scaling for the state vectors
        %
        % Can either be a scalar that will be used for every `f`-dimension, or a vector of the same
        % dimension as the system's core function `f` which will then be applied component-wise
        %
        % See the @ref scaling documentation for more details.
        %
        % @propclass{scaling}
        %
        % @type colvec
        StateScaling = 1;
        
        % The global sparsity pattern for the entire RHS
        SparsityPattern;
    end
    
    properties(SetAccess=protected, Transient)
        % The current parameter `\mu` for simulations, [] is none used.
        mu = [];
        
        % The current input function `u(t)` as function handle, [] if none used.
        u = [];
        
        % The current inputindex of the function `u(t)`
        inputidx = [];
    end
    
    properties(SetAccess=protected)
        NumStateDofs = [];
        NumAlgebraicDofs = 0;
        NumTotalDofs = [];
    end
    
    properties(SetAccess=private, Dependent)
        % The number of inputs available.
        InputCount;
        
        % The number of the system's parameters.
        ParamCount;
    end
    
    properties(SetAccess=private)
        % The Model this System is attached to.
        Model;
    end
        
    methods      
        function this = BaseFirstOrderSystem(model)
            % Creates a new base dynamical system class instance.
            %
            % Uses default output mapping: state variables are output.
            this.validateModel(model);
            this.Model = model;
            this.C = dscomponents.LinearOutputConv(1);
            
            % Register default properties
            this.registerProps('A','f','B','C','x0','Inputs','Params','MaxTimestep','StateScaling');
        end
        
        function rsys = getReducedSystemInstance(~, rmodel)
            % Creates a reduced system given the current system and the
            % reduced model.
            rsys = models.ReducedSystem(rmodel);
        end
        
        function setConfig(this, mu, inputidx)
            % Sets the dynamical system's configuration
            %
            % With configuration are meant the parameter `\mu` and input `u(t)` that effectively
            % results in a different dynamical system.
            %
            % Parameters:
            % mu: The current parameter `\mu\in\P`. Must be a column vector of the same length as
            % configured parameters. Leave empty if none are used.
            % inputidx: The index of the input function `u(t)` to use. Must be a valid index for one
            % of the cell elements of the Inputs property. Leave empty of none are used.
            %
            % See also: Inputs Params setParam addParam
            if isempty(mu) && ~isempty(this.Params)
                error('A model with parameters cannot be simulated without a parameter argument.');
            end
            if size(mu,1) ~= this.ParamCount
                error('The mu vector size mismatches the defined parameter number.');    
            end
            this.mu = mu;
            
            if isempty(inputidx) && ~isempty(this.Inputs)
                cprintf('red','Attention: Starting simulations without input, but the system has configured some.\n');
                %inputidx = 1;
            elseif inputidx > length(this.Inputs)
                error('Invalid input index "%d". There are %d possible inputs configured.',inputidx,length(this.Inputs));
            end
            this.u = @(t)[];
            this.inputidx = [];
            if ~isempty(inputidx)
                this.u = this.Inputs{inputidx};
                this.inputidx = inputidx;
            end
        end
        
        function prepareSimulation(this, mu, inputidx)
            this.setConfig(mu, inputidx);
            
            if isempty(this.NumStateDofs)
                this.NumStateDofs = size(this.x0.evaluate(mu),1);
                fprintf(2,'NumStateDofs not set. Guess from initial value: %d\n',this.NumStateDofs);
                this.updateDimensions;
            end
            
            % Forward preparation call to nonlinearity, if present
            if ~isempty(this.A)
                this.A.prepareSimulation(mu);
            end
            
            % Forward preparation call to nonlinearity, if present
            if ~isempty(this.f)
                this.f.prepareSimulation(mu);
            end
            
            % Forward preparation call to nonlinearity, if present
            if ~isempty(this.B)
                this.B.prepareSimulation(mu);
            end
        end
        
        function dx = ODEFun(this,t,x)
            % Debug variant for single evaluation. Commented in function
            % above.
            dx = zeros(this.NumTotalDofs,1);
            if ~isempty(this.A)
                dx = dx + this.A.evaluate(x, t);
            end
            if ~isempty(this.f)
                dx = dx + this.f.evaluate(x, t);
            end
            if ~isempty(this.B) && ~isempty(this.inputidx)
                dx = dx + this.B.evaluate(t, this.mu)*this.u(t);
            end
            if ~isempty(this.g)
                dx = [dx; this.g.evaluate(x,t)];
            end
        end
        
        function y = computeOutput(this, x, mu)
            % Computes the output `y(t) = C(t,\mu)Sx(t)` from a given state
            % result vector `x(t)`, using the system's time and current mu (if given).
            %
            % The matrix `S` represents possibly set state space scaling, and `C(t,\mu)` is an
            % output conversion. Identity values are assumed if the corresponding components
            % are not set/given.
            %
            % Parameters:
            % x: The state variable vector at each time step per column @type matrix<double>
            % mu: The parameter to use. If not specified, the currently set
            % `\mu` value of this system is used. @type colvec<double>
            % @default this.mu
            %
            % Return values:
            % y: The output according to `y(t) = C(t,\mu)Sx(t)`
            %
            % NOTE: This is also called for reduced simulations within ReducedModel.simulate.
            % However, reduced models do not employ state scaling anymore as it has been
            % included in the x0, C components at build-time for the reduced model,
            % respectively. Consequently, the StateScaling property of ReducedSystems is 1.
            %
            % See models.ReducedSystem.setReducedModel
            
            if nargin < 3
                mu = this.mu;
            end
            
            % Re-scale state variable
            if ~isequal(this.StateScaling,1)
                if isscalar(this.StateScaling)
                    x = this.StateScaling*x;
                else
                    x = bsxfun(@times,x,this.StateScaling);
                end
            end
            y = x;
            if ~isempty(this.C)
                if this.C.TimeDependent
                    % Evaluate the output conversion at each time t
                    % Figure out resulting size of C*x evaluation
                    t = this.Model.scaledTimes;
                    hlp = this.C.evaluate(t(1),mu)*x(:,1);
                    y = zeros(size(hlp,1),length(t));
                    y(:,1) = hlp;
                    for idx=2:length(t)
                        y(:,idx) = this.C.evaluate(t(idx),mu)*x(:,idx);
                    end
                else
                    % otherwise it's a constant matrix so multiplication
                    % can be preformed much faster.
                    y = this.C.evaluate([],mu)*x;
                end
            end
        end
        
        function updateSparsityPattern(this)
            td = this.NumTotalDofs;
            sd = this.NumStateDofs;
            JP = logical(sparse(td,td));
            if ~isempty(this.A) && ~isempty(this.A.JSparsityPattern)
                JP(1:sd,1:sd) = JP(1:sd,1:sd) | this.A.JSparsityPattern;
            end
            if ~isempty(this.f) && ~isempty(this.f.JSparsityPattern)
                JP(1:sd,:) = JP(1:sd,:) | this.f.JSparsityPattern;
            end
            if ~isempty(this.g)
                JP(sd+1:end,:) = this.g.JSparsityPattern;
            end
            this.SparsityPattern = JP;
        end
        
        function J = getJacobian(this, t, xc)
            % Computes the global jacobian of the current system.
            td = this.NumTotalDofs;
            sd = this.NumStateDofs;
            x = xc(1:sd);
            J = sparse(td,td);
            if isempty(this.g)
                if ~isempty(this.A)
                    J = J + this.A.getStateJacobian(x, t);
                end
                if ~isempty(this.f)
                    J = J + this.f.getStateJacobian(xc, t);
                end
            else
                if ~isempty(this.A)
                    J(1:sd,1:sd) = J(1:sd,1:sd) + this.A.getStateJacobian(x, t);
                elseif ~isempty(this.f)
                    J(1:sd,:) = J(1:sd,:) + this.f.getStateJacobian(xc, t);
                end
                J(sd+1:end,:) = this.g.getStateJacobian(xc,t);
            end
        end
        
        function x0 = getX0(this, mu)
            % Gets the initial state variable at `t=0`.
            %
            % This is exported into an extra function as it gets overridden
            % in the ReducedModel subclass, where ErrorEstimators possibly
            % change the x0 dimension.
            %
            % Parameters:
            % mu: The parameter `\mu` to evaluate `x_0(\mu)`. Use [] for
            % none. @typerowvec<double> @default []
            %
            % Return values:
            % x0: the initial state for the parameter `\mu` @type matrix<double>
            ss = this.StateScaling;
            if (size(mu,2) > 1) && ~isscalar(ss)
                ss = repmat(ss,1,size(mu,2));
            end
            x0 = this.x0.evaluate(mu) ./ ss;
            
            % Append DoFs for algebraic conditions below (if set)
            if this.NumAlgebraicDofs > 0
                x0 = [x0; this.getAlgebraicDofsInitialConditions];
            end
        end
        
        function M = getMassMatrix(this)
            % For first order systems, only algebraic constraints need to
            % be catered for.
            M = this.M;
            if isempty(M)
                M = dscomponents.ConstMassMatrix(speye(this.NumStateDofs));
            else
                if this.NumAlgebraicDofs > 0
                    error('First order with AlgebraicDofs not yet implemented.');
                end
            end
        end
        
        function p = addParam(this, name, default, varargin)
            % Adds a parameter with the given values to the parameter
            % collection of the current dynamical system.
            %
            % Use in subclass constructors to easily define desired default
            % parameters for a specific dynamical system.
            %
            % Parameters:
            % name: The name of the Parameter @type char
            % range: The range of the Parameter. Can be either a scalar or
            % a 1x2 double vector.
            % desired: The desired number of samples for that parameter.
            % @type integer @default 1
            % spacing: The intended sample spacing over the range @type
            % char @default 'lin'
            %
            % Return values:
            % p: The new ModelParam instance. @type ModelParam
            % 
            % See also: ModelParam setParam
            if ~isempty(this.Params) && ~isempty(find(strcmpi(name,{this.Params(:).Name}),1))
                error('Parameter with name %s already exists. Use setParam instead.',name);
            end
            p = data.ModelParam(name, default, varargin{:});
            this.Params(end+1) = p;
        end
        
        function pidx = getParamIndexFromName(this, paramname)
            % Gets the index within the parameter vector `\mu` for a given parameter name.
            %
            % Return values:
            % pidx: The parameter index for the given name, [] if none found.
            if ~ischar(paramname)
                error('Parameter paramname must be a char array');
            end
            pidx = find(strcmpi(paramname,{this.Params(:).Name}),1);
        end
        
        function pt = getParamInfo(this, mu)
            if nargin < 2
                mu = this.Model.DefaultMu;
            end
            %str = '';
            pt = PrintTable;
            pt.HasRowHeader = true;
            pt.HasHeader = true;
            pt.addRow('Parameter name','Value','Index');
            for i=1:this.ParamCount
               % str = sprintf('%s%d - %s: %f\n',str,i,this.Params(i).Name,mu(i));
                pt.addRow(this.Params(i).Name,mu(i),i);
            end
            if nargout < 1
                pt.print;
            end
        end
        
        function plotInputs(this, pm)
            if nargin < 2
                pm = PlotManager;
                pm.LeaveOpen = true;
            end
            for i=1:this.InputCount
                ui = this.Inputs{i};
                h = pm.nextPlot(sprintf('u%d',i),sprintf('u_%d = %s',i,func2str(ui)),...
                    't',sprintf('u_%d',i));
                plot(h,this.Model.Times,ui(this.Model.scaledTimes),'LineWidth',2);
            end
            pm.done;
        end
    end
    
    methods(Access=protected)
        
        function updateDimensions(this)
            this.NumTotalDofs = this.NumStateDofs + this.NumAlgebraicDofs;
        end
        
        function ad_ic = getAlgebraicDofsInitialConditions(this)
            % The default is to return all zeros.
            ad_ic = zeros(this.NumAlgebraicDofs,1);
        end
        
        function validateModel(this, model)%#ok
            % Validates if the model to be set is a valid BaseModel at
            % least.
            % Extracting this function out of the setter enables subclasses
            % to further restrict the objects that may be passed, as is
            % being done in models.ReducedSystem, for example.
            if ~isa(model, 'models.BaseModel')
                error('The Model property has to be a child of models.BaseModel');
            end
        end
    end
        
    %% Getter & Setter
    methods
        function set.f(this,value)
            if ~isempty(value) && ~isa(value, 'dscomponents.ACoreFun')
                error('The property "f" has to be a class implementing dscomponents.ACoreFun');
            end
%             if (isequal(this,value))
%                 warning('KerMor:selfReference','Careful with self-referencing classes. See BaseFullModel class documentation for details.');
%             end
            this.f = value;
        end
        
        function set.B(this,value)
            if ~isempty(value) && ~isa(value, 'dscomponents.AInputConv')
                error('The property "B" has to be a class implementing dscomponents.AInputConv');
            end
            this.B = value;
        end
        
        function set.A(this,value)
            if ~isempty(value) && ~(isa(value, 'dscomponents.LinearCoreFun') || isa(value, 'dscomponents.AffLinCoreFun'))
                error('The property "A" has to be a LinearCoreFun or AffLinCoreFun.');
            end
            this.A = value;
        end
        
        function set.C(this,value)
            if isempty(value)
                value = dscomponents.LinearOutputConv(1);
                warning('An output conversion must always exist. Choosing dscomponents.LinearOutputConv(1) for simple forwarding.');
            elseif ~isa(value, 'dscomponents.AOutputConv')
                error('The property "C" has to be a class implementing dscomponents.AOutputConv');
            end
            this.C = value;
        end
        
        function set.Inputs(this, value)
            if ~iscell(value)
                error('Property "Inputs" must be a cell array.');
            end
            for n=1:numel(value)
                if ~isempty(value{n})
                    if ~isa(value{n},'function_handle')
                        error('Each "Inputs" cell must contain a function handle.');
%                     elseif nargin(value{n}) ~= 1
%                        error('Each "Inputs" function must take exactly one (=time) parameter.');
                    end
                end
            end
            this.Inputs = value;
        end
                
        function set.Params(this, value)
            if ~isa(value,'data.ModelParam')
                error('Params property must be a ModelParam array.');
            end
            this.Params = value;
        end
        
        function set.x0(this, value)
            this.checkType(value,'dscomponents.AInitialValue');
%             if ~isa(value,'function_handle')
%                 error('x0 must be a function handle.');
%             elseif nargin(value) ~= 1
%                 error('x0 must take exactly one argument.');
%             end
            
            this.x0 = value;
        end
        
        function value = get.ParamCount(this)
            value = length(this.Params);
        end
        
        function value = get.InputCount(this)
            value = length(this.Inputs);
        end
        
        function set.MaxTimestep(this, value)
            if (~isscalar(value) || value <= 0) && ~isempty(value)
                error('Value must be a positive real scalar if not empty.');
            end
            this.MaxTimestep = value;
        end
        
        function set.StateScaling(this, value)
            if isempty(value)
                error('StateScaling must not be empty. (Set to 1 if not needed).');
            elseif ~isvector(value)
                error('StateScaling must be a vector.');
            end
            this.StateScaling = value(:);
        end
    end
        
end

