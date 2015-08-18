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
            this.x0deriv = fullsys.x0deriv.project(rm.V,rm.W);
            % Project damping matrix
            if ~isempty(fullsys.D)
                this.D = fullsys.D.project(rm.V,rm.W);
            end
        end
        
        function J = getJacobian(this, t, x_xdot_c)
            % Computes the global jacobian of the current RHS system.
            this.nJevals = this.nJevals + 1;
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
        
        function dx = ODEFun(this,t,x)
            est = this.Model.ErrorEstimator;
            haveest = ~isempty(est) && est.Enabled;
            if haveest
                xall = x;
                x = x(1:end-est.ExtraODEDims,:);    
            end
            
            dx = ODEFun@models.BaseSecondOrderSystem(this,t,x);
            
            if haveest
                dx = [dx; est.evalODEPart(xall, t, this.u(t))];
            end
        end
        
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
            
            % Insert reduced derivative initial values (xdot0) between
            z_zdot_c0(num_z_dof+(1:num_zdot_dof)) = this.x0deriv.evaluate(mu);
            
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
        
        function val = getDerivativeDirichletValues(this, t)
            fsys = this.Model.FullModel.System;
            val = fsys.getDerivativeDirichletValues(t);
        end
        
    end
    
    methods(Access=protected)
        
        function R = compileReconstructionMatrix(this, V)
            if any(this.DerivativeDirichletPosInStateDofs)
                sd = this.NumStateDofs;
                dd = this.NumDerivativeDofs;
                numdd = sd-dd;
                fsys = this.FullSystem;
                Vstate = zeros(fsys.NumStateDofs,this.NumStateDofs);
                Vstate(fsys.DerivativeDirichletPosInStateDofs,1:numdd) = eye(numdd);
                Vstate(~fsys.DerivativeDirichletPosInStateDofs,numdd+1:end) = V;
            else
                Vstate = V;
            end
            R = blkdiag(Vstate,V,eye(this.NumAlgebraicDofs));
        end
        
        function x0proj = projectx0(this,x0,SV,SW)
            if any(this.DerivativeDirichletPosInStateDofs)
                fsys = this.FullSystem;
                SWD = zeros(fsys.NumStateDofs,this.NumStateDofs);
                numdd = fsys.NumStateDofs - fsys.NumDerivativeDofs;
                SWD(fsys.DerivativeDirichletPosInStateDofs,1:numdd) = eye(numdd);
                SWD(~fsys.DerivativeDirichletPosInStateDofs,numdd+1:end) = SW;
                x0proj = x0.project(SV,SWD);
            else
                x0proj = projectx0@models.ReducedSystem(this,x0,SV,SW);
            end
        end
        
        function updateDimensions(this)
            updateDimensions@models.ReducedSystem(this);
            this.NumDerivativeDofs = this.NumStateDofs;
            
            % Take care of extra derivative dirichlet values - they are not
            % reduced but kept as full variables in the reduced system (not
            % worth the effort)
            fsys = this.Model.FullModel.System;
            numdd = fsys.NumStateDofs - fsys.NumDerivativeDofs;
            if numdd > 0
                this.NumStateDofs = this.NumStateDofs + numdd;
                pos = false(this.NumStateDofs,1);
                pos(1:numdd) = true;
                this.DerivativeDirichletPosInStateDofs = pos;
            end
            
            this.NumTotalDofs = this.NumStateDofs + ...
                this.NumDerivativeDofs + this.NumAlgebraicDofs;
        end
        
        function validateModel(this, model)%#ok
            validateModel@models.ReducedSystem(this, model);
        end
    end
    
end

