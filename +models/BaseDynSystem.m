classdef BaseDynSystem < handle
    %BASEDYNSYSTEM Base class for all KerMor dynamical systems.
    %
    % To setup custom dynamical systems, inherit from this class.
    %
    % Programming convention:
    % Any system-specific properties should be added in the subclass, even
    % though sometimes one might tend to put settings into the core
    % function f (e.g. for space-discretized systems and thus all
    % space-discretization relevant settings). In order to have all
    % settings combined at one class one should rather pass the custom
    % system's instance to the ACoreFun implementing classes' constructor
    % and read settings from that instead.
    % For an example see the models.pcd.PCDSystem and models.pcd.CoreFun.
    %
    % @author Daniel Wirtz @date 17.03.2010
    
    properties
        % The core f function from the dynamical system.
        %
        % @type dscomponents.ACoreFun
        f;
        
        % The input conversion
        B;
        
        % The output conversion
        % Defaults to an LinearOutputConv instance using a 1-matrix, which
        % just forwards the state variables and supports projection.
        %
        % @sa LinearOutputConv::project(V)
        C;
        
        % See also: dscomponents.LinearOutputConv/project
        
        % The system's possible input functions.
        % A cell array of function handles, each taking a time argument t.
        Inputs = {};
        
        % The parameters usable for the dynamical system.
        Params = models.ModelParam.empty;
        
        % The initial value function.
        % A function handle taking the parameter argument mu.
        % The argument mu may be empty if no parameters are used within the
        % dynamic system.
        x0;
        
        % The maximum timestep allowed for any ODE solvers.
        %
        % This might be necessary if the Core function encapsulates a
        % spatial discretization and thus CFL conditions apply, for example.
        %
        % @note If any time scaling is used within the model, this value
        % must correspond to the ''scaled'' maximum timestep.
        %
        % See also: BaseModel BaseModel.tau BaseModel.dtscaled
        MaxTimestep = [];
    end
    
    properties(Dependent)
        % The number of inputs available.
        InputCount;
        
        % The number of the system's parameters.
        ParamCount;
    end
    
    methods
        
        function this = BaseDynSystem
            % Creates a new base dynamical system class instance.
            %
            % Uses default output mapping: state variables are output.
            % Initial conditions are set to zero.
            %
            % Parameters:
            % model: The parent model that 
            this.C = dscomponents.LinearOutputConv(1);
            this.x0 = @(mu)0;
        end
        
        function odefun = getODEFun(this, mu, inputidx)
            % Generates the ODE function for a specific parameter mu and
            % given input index.
            %
            % Depending on the system mu may be [] and is forwarded to the
            % evaluate function of the system's components.
            % Analogue the inputidx may be empty to force skipping the
            % system's input; however, if no inputconv component is set
            % this is ignored in any case.
                        
            % System without inputs
            if this.InputCount == 0 || isempty(this.B)
                odefun = @(t,x)(this.f.evaluate(x,t,mu));
            else
                % generates the ode function for given parameter and input function
                u = this.Inputs{inputidx};
                odefun = @(t,x)(this.f.evaluate(x,t,mu) + this.B.evaluate(t,mu)*u(t));
            end
        end
        
        function mu = getRandomParam(this)
            % Gets a random parameter sample from the system's parameter
            % domain P
            if this.ParamCount > 0
                pmin = [this.Params(:).MinVal]';
                pmax = [this.Params(:).MaxVal]';
                mu = rand(this.ParamCount,1) .* (pmax-pmin) + pmin;
            else
                mu = [];
            end
        end
        
    end
    
    %% Getter & Setter
    methods
        
        function addParam(this, name, range, desired)
            % Adds a parameter with the given values to the parameter
            % collection of the current dynamical system.
            %
            % Use in subclass constructors to easily define desired default
            % parameters for a specific dynamical system.
            %
            % Parameters:
            % name: The name of the Parameter
            % range: The range of the Parameter. Can be either a scalar or
            % a 1x2 double vector.
            % desired: The desired number of samples for that parameter.
            % Defaults to 1.
            %
            % See also: ModelParam setParam
            
            if nargin < 4
                desired = 1;
            end
            if ~isempty(this.Params) && ~isempty(find(strcmpi(name,{this.Params(:).Name}),1))
                error('Parameter with name %s already exists. Use setParam instead.',name);
            end
            this.Params(end+1) = models.ModelParam(name, range, desired);
        end
        
        function setParam(this, name, range, desired)
            % Sets values for a parameter with the name "name".
            % If no such parameter exists a new one is created using the
            % specified name and values.
            %
            % Parameters:
            % name: The name of the Parameter
            % range: The range of the Parameter. Can be either a scalar or
            % a 1x2 double vector.
            % desired: The desired number of samples for that parameter.
            % Optional; defaults to 1.
            %
            % See also: addParam
            
            if nargin < 4
                desired = 1;
            end
            if ~isempty(this.Params)
                pidx = find(strcmpi(name,{this.Params(:).Name}),1);
                if ~isempty(pidx)
                    this.Params(pidx).Range = range;
                    if ~isempty(desired)
                        this.Params(pidx).Desired = desired;
                    end
                    return;
                end
            end
            this.addParam(name, range, desired);
        end
        
        function set.f(this,value)
            if ~isempty(value) && ~isa(value, 'dscomponents.ACoreFun')
                error('The property "f" has to be a class implementing dscomponents.ACoreFun');
            end
            this.f = value;
        end
        
        function set.B(this,value)
            if ~isempty(value) && ~isa(value, 'dscomponents.IInputConv')
                error('The property "B" has to be a class implementing dscomponents.IInputConv');
            end
            this.B = value;
        end
        
        function set.C(this,value)
            if ~isempty(value) && ~isa(value, 'dscomponents.AOutputConv')
                error('The property "C" has to be a class implementing dscomponents.AOutputConv');
            end
            this.C = value;
        end
        
        function set.x0(this,value)
            if ~isa(value,'function_handle')
                error('Argument funPtr must be a function handle.');
            elseif nargin(value) ~= 1
                error('funPtr nargin must equal one (= mu).');
            end
            this.x0 = value;
        end
        
        function set.Inputs(this, value)
            if ~iscell(value)
                error('Property "Inputs" must be a cell array.');
            end
            for n=1:length(value)
                if ~isa(value{n},'function_handle')
                    error('Each "Inputs" cell must contain a function handle.');
                elseif nargin(value{n}) ~= 1
                    error('Each "Inputs" function must take exactly one (=time) parameter.');
                end
            end
            this.Inputs = value;
        end
        
        function set.Params(this, value)
            if ~isa(value,'models.ModelParam');
                error('Params property must be a ModelParam array.');
            end
            this.Params = value;
        end
        
        function value = get.ParamCount(this)
            value = length(this.Params);
        end
        
        function value = get.InputCount(this)
            value = length(this.Inputs);
        end
    end
    
end
