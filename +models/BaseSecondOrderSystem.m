classdef BaseSecondOrderSystem < models.BaseFirstOrderSystem
    % Base class for all KerMor second-order dynamical systems.
    %
    % To setup custom second-order dynamical systems, inherit from this class.
    %
    % @author Daniel Wirtz @date 2010-03-17
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetObservable)
        % The damping matrix of the second order system.
        %
        % @propclass{important}
        %
        % @type dscomponents.AffLinCoreFun,dscomponents.LinearCoreFun
        % @default []
        D;
    end
    
    properties(SetAccess=protected)
        NumSecondOrderDoFs;
    end
        
    methods      
        function this = BaseSecondOrderSystem(model)
            % Creates a new base dynamical system class instance.
            this = this@models.BaseFirstOrderSystem(model);
            
            % Register default properties
            this.registerProps('D');
        end
        
%         function createComponents(this)
%             createComponents@models.BaseFirstOrderSystem(this);
%             D_ = this.getD;
%             if ~isempty(D_)
%                 this.D = D_;
%             end
%         end
        
        function rsys = buildReducedSystem(this, rmodel)%#ok
            % Creates a reduced system given the current system and the
            % reduced model.
            %
            % Override this method in subclasses (with call to parent) if
            % custom behaviour of your system is required prior/posterior
            % to reduced model computation.
            error('TODO');
            rsys = models.ReducedSystem(rmodel);
        end
        
        function prepareSimulation(this, mu, inputidx)
            prepareSimulation@models.BaseFirstOrderSystem(this, mu, inputidx);
            
            % Forward preparation call to D, if present
            if ~isempty(this.D)
                this.D.prepareSimulation(mu);
            end
        end
    
        function odefun = getODEFun(this)
            odefun = @this.ODEFun;
            return;
%             % Determine correct ODE function (A,f,B combination)
%             str = {};
%             if ~isempty(this.A)
%                 str{end+1} = 'this.A.evaluate(x, t)';
%             end
%             if ~isempty(this.f)
%                 str{end+1} = 'this.f.evaluate(x, t)';
%             end
%             if ~isempty(this.B) && ~isempty(this.inputidx)
%                 str{end+1} = 'this.B.evaluate(t, this.mu)*this.u(t)';
%             end
%             odefunstr = Utils.implode(str,' + ');
%             if ~isempty(this.g)
%                 odefunstr = ['[' odefunstr '; this.g.evaluate(x,t)]'];
%             end
%             odefun = eval(['@(t,x)' odefunstr]);
        end
        
        function dxyc = ODEFun(this, t, xyc)
            % The state space vector is composed of
            % x: Original state space of second order model
            % y: Substituted variable x'=y for first order transformation
            % c: Algebraic constraint variables
            dx = zeros(size(xyc));
            
            odof = this.NumSecondOrderDoFs;
            x = xyc(odof+1:2*odof);
            y = xyc(1:odof);
            
            % Set x'=y
            dx(1:odof) = x;
            
            % set y' = A(x)+D(y)+f(x,y)+Bu
            dy = zeros(odof,1);
            if ~isempty(this.A)
                dy = dy + this.A.evaluate(x, t);
            end
            if ~isempty(this.D)
                dy = dy + this.D.evaluate(y, t);
            end
            if ~isempty(this.f)
                dy = dy + this.f.evaluate(xyc, t);
            end
            if ~isempty(this.B) && ~isempty(this.inputidx)
                B = this.B.evaluate(t, this.mu)*this.u(t);
                dy = dy + B;
            end
            dxyc = [dx;dy];
            % Append algebraic constraint evaluation
            if ~isempty(this.g)
                dxyc = [dxyc; this.g.evaluate(xyc,t)];
            end
        end
        
        function [JacFun, JPattern] = getJacobianInfo(this)
            d = this.NumSecondOrderDoFs;
            ad = 0;
            if ~isempty(this.g)
                ad = this.g.fDim;
                %length(this.AlgebraicConditionDoF);
            end 
            JTrans = [sparse(d,d) speye(d) sparse(d,ad)];
            JacFunStr = '[JTrans; ';
            JPattern = JTrans;
            if ~isempty(this.A) && ~isempty(this.f)
                JacFunStr = [JacFunStr 'this.A.getStateJacobian(x, t) + this.f.getStateJacobian(x, t);'];
                if ~isempty(this.A.JSparsityPattern) && ~isempty(this.f.JSparsityPattern)
                    [i,j] = find(this.A.JSparsityPattern + this.f.JSparsityPattern);
                    JPattern = [JPattern; sparse(i,j,ones(length(i),1),...
                        this.A.fDim,this.A.xDim)];
                end
            elseif ~isempty(this.A)
                JacFunStr = [JacFunStr 'this.A.getStateJacobian(x, t);'];
                if  ~isempty(this.A.JSparsityPattern)
                    JPattern = [JPattern; this.A.JSparsityPattern];
                end
            elseif ~isempty(this.f)
                JacFunStr = [JacFunStr 'this.f.getStateJacobian(x, t);'];
                if  ~isempty(this.f.JSparsityPattern)
                    JPattern = [JPattern; this.f.JSparsityPattern];
                end
            end
            % Algebraic dofs
            if ~isempty(this.g)
                JPattern = [JPattern; this.g.JSparsityPattern];
                JacFunStr = [JacFunStr 'this.g.getStateJacobian(x,t)'];
            end
            spy(JPattern)
            JacFunStr
            JacFun = eval(['@(t,x)' JacFunStr ']']);
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
    
    methods(Sealed)
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
    end
    
    methods(Access=protected)
        
        function x0 = getX0(this, mu)
            % Gets the initial state variable at `t=0`.
            x0 = getX0@models.BaseFirstOrderSystem(this, mu);
            odof = this.NumSecondOrderDoFs - length(this.AlgebraicConditionDoF);
            x0 = [zeros(odof,1); x0];
            
            % TODO:
            % Remove dirichlet values
            %x0(this.idx_uv_bc_glob) = [];
        end
        
        function gx = computeAlgebraicConditions(this, x, t)
            % Callback to compute algebraic conditions, if any.
            gx = [];
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
end

