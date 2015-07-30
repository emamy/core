classdef ReducedSecondOrderSystem < models.ReducedSystem & models.BaseSecondOrderSystem
    %REDUCEDSECONDORDERSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function this = ReducedSecondOrderSystem(rmodel)
            % Creates a new base dynamical system class instance.
            this = this@models.BaseSecondOrderSystem(rmodel);
            this = this@models.ReducedSystem(rmodel);
        end
        
        function build(this)
            build@models.ReducedSystem(this);
            rm = this.Model;
            fullsys = rm.FullModel.System; 
            % Additional steps
            % Project damping matrix
            if ~isempty(fullsys.D)
                this.D = fullsys.D.project(rm.V,rm.W);
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
            I = speye(sd);
            I(:,this.DerivativeDirichletPosInStateDofs) = [];
            J(xpos,:) = [sparse(sd,sd) I sparse(sd,ad)];
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
        
%         function odefun = getODEFun(this)
%             % Determine correct ODE function (A,f,B combination)
%             xarg = 'x';
%             est = this.Model.ErrorEstimator;
%             haveest = ~isempty(est) && est.Enabled;
%             if haveest
%                 xarg = 'x(1:end-est.ExtraODEDims,:)';
%             end
%             
%             str = {};
%             if ~isempty(this.A)
%                 str{end+1} = sprintf('this.A.evaluate(%s, t)',xarg);
%             end
%             if ~isempty(this.f)
%                 str{end+1} = sprintf('this.f.evaluate(%s, t)',xarg);
%             end
%             if ~isempty(this.B) && ~isempty(this.inputidx)
%                 str{end+1} = 'this.B.evaluate(t, this.mu)*this.u(t)';
%             end
%             funstr = Utils.implode(str,' + ');
%             if haveest
%                 funstr = ['[' funstr '; est.evalODEPart(x, t, this.u(t))]'];
%             end
%             odefun = eval(['@(t,x)' funstr]);
%         end
        
        function z_zdot_c0 = getX0(this, mu)
            % Gets the initial value of `x_0(\mu)`.
            %
            % If the estimator is enabled, x0 is extended by the e0
            % components of the error estimator.
            
            z_c0 = getX0@models.BaseFirstOrderSystem(this, mu);
            num_z_dof = this.NumStateDofs;
            num_zdot_dof = this.NumDerivativeDofs;
            z_zdot_c0 = zeros(this.NumTotalDofs,1);
            
            % Insert state initial values
            z_zdot_c0(1:num_z_dof) = z_c0(1:num_z_dof);
            
            % Insert zero derivative initial values (xdot0) between
            z_zdot_c0(num_z_dof+(1:num_zdot_dof)) = zeros(num_zdot_dof,1);
            
            % Insert alg cond initial values
            z_zdot_c0(num_z_dof+num_zdot_dof+1:end) = z_c0(num_z_dof+1:end);
            
            m = this.Model;
            if ~isempty(m.ErrorEstimator) && m.ErrorEstimator.Enabled
                z_zdot_c0 = [z_zdot_c0; m.ErrorEstimator.getE0(mu)];
            end
        end
        
        function updateSparsityPattern(this)
            updateSparsityPattern@models.ReducedSystem(this);
        end
        
    end
    
    methods(Access=protected)
        
        function R = compileReconstructionMatrix(this, V)
            R = blkdiag(V,V,eye(this.NumAlgebraicDofs));
        end
        
        function updateDimensions(this)
            updateDimensions@models.ReducedSystem(this);
            this.NumDerivativeDofs = this.NumStateDofs;
            this.NumTotalDofs = 2*this.NumStateDofs + this.NumAlgebraicDofs;
        end
        
        function val = getDerivativeDirichletValues(this, t)
            val = getDerivativeDirichletValues@models.BaseSecondOrderSystem(this, t);
        end
        
        function validateModel(this, model)%#ok
            validateModel@models.ReducedSystem(this, model);
        end
    end
    
end

