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
        NumDerivativeDofs;
        
        % A logical column vector containing true at the locations of
        % explicit derivative dirichlet conditions.
        % The vector must have the same size as NumStateDofs.
        %
        % See also: getDerivativeDirichletValues
        DerivativeDirichletPosInStateDofs;
    end
        
    properties(Access=protected)
        % Pre-Computed pattern for dx-part
        JDX;
    end
    
    methods      
        function this = BaseSecondOrderSystem(model)
            % Creates a new base dynamical system class instance.
            this = this@models.BaseFirstOrderSystem(model);
            
            % Register default properties
            this.registerProps('D');
        end
        
        function rsys = getReducedSystemInstance(~, rmodel)
            % Creates a reduced system given the current system and the
            % reduced model.
            %
            % Override this method in subclasses (with call to parent) if
            % custom behaviour of your system is required prior/posterior
            % to reduced model computation.
            rsys = models.ReducedSecondOrderSystem(rmodel);
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
        
        function dx_xdot_c = ODEFun(this, t, x_xdot_c)
            % The state space vector is composed of
            % x: Original state space of second order model
            % xdot: Substituted variable x'=y for first order transformation
            % c: Algebraic constraint variables
            dx_xdot_c = zeros(size(x_xdot_c));
            
            num_x_dof = this.NumStateDofs;
            num_xdot_dof = this.NumDerivativeDofs;
            x = x_xdot_c(1:num_x_dof);
            xdot = x_xdot_c(num_x_dof+(1:num_xdot_dof));
            
            %% Set x'=xdot
            % Check for derivative dirichlet conditions - insert them here
            % as direct values
            if num_x_dof ~= num_xdot_dof % Quickest check
                dx = zeros(num_x_dof,1);
                % Execute callback - subclasses must provide the (possibly
                % time-dependent) valuess
                val = this.getDerivativeDirichletValues(t);
                pos = this.DerivativeDirichletPosInStateDofs;
                % Dirichlet values
                dx(pos) = val;
                % Dof values
                dx(~pos) = xdot;
                % Set in global dof vector
                dx_xdot_c(1:num_x_dof) = dx;
            else
                % Direct forwarding - no derivative dirichlet conditions
                dx_xdot_c(1:num_x_dof) = xdot;
            end
            
            %% Set x''=xdot' = A(x)+D(y)+f(x,y)+Bu
            dy = zeros(num_xdot_dof,1);
            if ~isempty(this.A)
                dy = dy + this.A.evaluate(x, t);
            end
            if ~isempty(this.D)
                dy = dy + this.D.evaluate(xdot, t);
            end
            if ~isempty(this.f)
                dy = dy + this.f.evaluate(x_xdot_c, t);
            end
            if ~isempty(this.B) && ~isempty(this.inputidx)
                B = this.B.evaluate(t, this.mu)*this.u(t);
                dy = dy + B;
            end
            dx_xdot_c(num_x_dof+(1:num_xdot_dof)) = dy;
            % Append algebraic constraint evaluation
            if ~isempty(this.g)
                dx_xdot_c(num_x_dof+num_xdot_dof+1:end) = this.g.evaluate(x_xdot_c,t);
            end
        end
        
        function J = getJacobian(this, t, x_xdot_c)
            % Computes the global jacobian of the current RHS system.
            td = this.NumTotalDofs;
            sd = this.NumStateDofs;
            dd = this.NumDerivativeDofs;
            ad = this.NumAlgebraicDofs;
            xdotpos = sd+(1:dd);
            xpos = 1:sd;
            cpos = sd+dd+(1:ad);
            x = x_xdot_c(xpos);
            xdot = x_xdot_c(xdotpos);
            J = sparse(td,td);
            J(xpos,:) = this.JDX;
            if ~isempty(this.A)
                J(xdotpos,xpos) = J(xdotpos,xpos) + this.A.getStateJacobian(x, t);
            end
            if ~isempty(this.D)
                J(xdotpos,xdotpos) = J(xdotpos,xdotpos) + this.D.getStateJacobian(xdot, t);
            end
            if ~isempty(this.f)
                J(xdotpos,:) = J(xdotpos,:) + this.f.getStateJacobian(x_xdot_c, t);
            end
            if ~isempty(this.g)
                J(cpos,:) = this.g.getStateJacobian(x_xdot_c,t);
            end
        end
        
        function updateSparsityPattern(this)
            % The state space vector (#NumTotalDofs) is composed of
            % x: Original state space of second order model, #NumStateDofs
            % xdot: Substituted variable x'=y for first order
            % transformation #NumDerivativeDoFs
            % c: Algebraic constraint variables #NumAlgebraicDofs
            td = this.NumTotalDofs;
            sd = this.NumStateDofs;
            dd = this.NumDerivativeDofs;
            ad = this.NumAlgebraicDofs;
            xdotpos = sd+(1:dd);
            JP = logical(sparse(td,td));
            I = speye(sd);
            I(:,this.DerivativeDirichletPosInStateDofs) = [];
            this.JDX = [sparse(sd,sd) I sparse(sd,ad)];
            JP(1:sd,:) = this.JDX;
            if ~isempty(this.A) && ~isempty(this.A.JSparsityPattern)
                JP(xdotpos,1:sd) = JP(xdotpos,1:sd) | this.A.JSparsityPattern;
            end
            if ~isempty(this.D) && ~isempty(this.D.JSparsityPattern)
                JP(xdotpos,xdotpos) = JP(xdotpos,xdotpos) | this.D.JSparsityPattern;
            end
            if ~isempty(this.f) && ~isempty(this.f.JSparsityPattern)
                JP(xdotpos,:) = JP(xdotpos,:) | this.f.JSparsityPattern;
            end
            if ~isempty(this.g)
                JP(sd+dd+1:end,:) = this.g.JSparsityPattern;
            end
            this.SparsityPattern = JP;
        end
        
        function x_xdot_c0 = getX0(this, mu)
            % Compiles the global x0 vector of the global dof vector
            
            % Gets the initial state variable at `t=0`.
            x_c0 = getX0@models.BaseFirstOrderSystem(this, mu);
            num_x_dof = this.NumStateDofs;
            num_xdot_dof = this.NumDerivativeDofs;
            x_xdot_c0 = zeros(this.NumTotalDofs,1);
            
            % Insert state initial values
            x_xdot_c0(1:num_x_dof) = x_c0(1:num_x_dof);
            
            % Insert zero derivative initial values (xdot0) between
            x_xdot_c0(num_x_dof+(1:num_xdot_dof)) = zeros(num_xdot_dof,1);
            
            % Insert alg cond initial values
            x_xdot_c0(num_x_dof+num_xdot_dof+1:end) = x_c0(num_x_dof+1:end);
        end
        
        function M = getMassMatrix(this)
            M = this.M;
            if ~isa(M,'dscomponents.ConstMassMatrix')
                error('Non-constant mass matrices not yet implemented for second order systems');
            end
            sd = this.NumStateDofs;
            ad = this.NumAlgebraicDofs;
            M = blkdiag(speye(sd),this.M.M,sparse(ad,ad));
            % Tell the mass matrix which components are algebraic
            % constraint dofs - those wont be used determinig a suitable
            % initial slope in e.g. solvers.MLWrapper
            algdofs = sd+this.NumDerivativeDofs+(1:ad);
            M = dscomponents.ConstMassMatrix(M,algdofs);
        end
    end
    
    methods(Access=protected)
        
        function updateDimensions(this)
            this.NumTotalDofs = this.NumStateDofs ...
                + this.NumDerivativeDofs + this.NumAlgebraicDofs;
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
    
    methods(Access=protected, Abstract)
        % Computes the derivative dirichlet values dependent on the
        % current time.
        %
        % The logical field DerivativeDirichletPosInStateDofs indicates
        % the positions of those values within the StateDofs.
        %
        % The default behaviour is doing nothing (just return []); this
        % function will only be called if NumStateDofs and
        % NumDerivativeDofs differ.
        %
        % See also: ODEFun DerivativeDirichletPosInStateDofs
        val = getDerivativeDirichletValues(this, t);
    end
end

