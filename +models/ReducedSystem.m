classdef ReducedSystem < models.BaseFirstOrderSystem
%ReducedSystem: A KerMor reduced dynamical system.
%
% Any state space scaling is included into the reduced model upon creation, as scaling is
% possibly of the large systems' dimension. See @ref scaling for more information.
% 
% @author Daniel Wirtz @date 17.03.2010
%
% @change{0,6,dw,2012-05-29}
% - New direct methods for ODE handles (selection once per call more
% effective than if clauses every evaluation)
% - fixed cases with A component
% - fixed setting of A component upon setReducedModel
%
% @change{0,5,dw,2011-10-15} Fixed the computation of reduced systems. Previously, the scaling
% was not included in the reduced model if no subspace projection was used. Now, scaling is
% included if set and independently from any subspace projection.
%
% @change{0,5,dw,2011-08-05} Fixed the use of state space scaling in reduced simulations. Now
% the scaling and reduction is performed as described in @ref state_scaling.
%
% @change{0,4,dw,2011-05-13} Removed the old 'getODEFun' function and replaced it by a direct
% implementation 'ODEFun'. This enables to avoid nested function handles with in turn allow for
% a speedup of reduced simulations by almost a factor of 2.
%
% See also: models.BaseFirstOrderSystem
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing    

    properties(SetAccess=protected)
        % The reduced dof to full dof reconstruction matrix
        R;
    end
    
    properties(Access=private)
        plotPtr;
    end
    
    properties(Dependent)
        FullSystem;
    end
    
    methods
        function this = ReducedSystem(rmodel)
            % Creates a new ReducedSystem instance.
            %
            % Parameters:
            % rmodel: [Optional] The reduced model to create the reduced
            % system from.
            this = this@models.BaseFirstOrderSystem(rmodel);
        end
        
        function build(this)
            % Creates a reduced system from BaseFullModel child system.
            % As default, all system's components are copied as-is and any
            % reduction changes are performed in the
            % BaseFullModel.buildReducedModel method.
            %
            % Parameters:
            % rmodel: The instance of the reduced model that will
            % contain this instance as System property.
            %
            % Important:
            % The reduced model must have set all reduced data (i.e. V,W)
            % that will be used during this constructor; references to
            % variables/properties from the full model will result in large
            % reduced models as references to the full models properties
            % also save the full model with all data.
            %
            % See also: models.BaseFullModel#buildReducedModel
            disp('Start building reduced system...');
            
            rmodel = this.Model;
            fullmodel = rmodel.FullModel;
            fullsys = fullmodel.System;
                        
            V = rmodel.V;
            
            if ~isempty(fullmodel.SpaceReducer) && isempty(V)
                error(['Model has a SpaceReducer but no projection'...
                        'matrix V is given. Forgot to call offlineGenerations?']);
            end
            
            W = rmodel.W;
            % Use V if W is empty (Galerkin projection)
            if ~isempty(V) && isempty(W)
                W = V;
            end
            
            this.updateDimensions;
                        
            % The state scaling for the reduced system is one, as the scaling matrices have been
            % incorporated into the V,W matrices via V := SV and W = S^-1W inside the ReducedModel
            % setModel method.
            this.StateScaling = 1;
            SV = 1; SW = 1; % Default: no scaling
            s = fullmodel.System.StateScaling;
            if s ~= 1
                if isscalar(s)
                    dim = fullsys.NumStateDofs;
                    s(1:dim,1) = s;
                else
                    dim = length(s);
                end
                SV = spdiags(s,0,dim,dim);
                SW = spdiags(1./s,0,dim,dim);
            end
            % Incorporate the state scaling into the projection matrices.
            % They are used to project the x0 initial values and output C.
            if ~isempty(V)
                SV = SV * V;
                SW = SW * W;
            end
            
            % Clones the full system's basic (all but functions)
            % properties
            this.Params = fullsys.Params;
            this.Inputs = fullsys.Inputs;
            this.MaxTimestep = fullsys.MaxTimestep;
            
            % Copy component handles (it's NOT cloning them!)
            % This is the default as if no reduction methods are
            % applied (sensless but yet allowed) an identical system
            % will be the result.
            this.f = fullsys.f;
            this.B = fullsys.B;
            
            % SV ~= 1 means that projection or nontrivial scaling is used.
            if ~isequal(SV,1) 
                if ~isempty(fullsys.x0)  % otherwise error for static models
                    this.x0 = this.projectx0(fullsys.x0,SV,SW);
%                     this.x0 = fullsys.x0.project(SV,SW);
                end
                SR = this.compileReconstructionMatrix(SV);
                this.C = fullsys.C.project(SR,SW);
            else
                if ~isempty(fullsys.x0)  % otherwise error for static models
                    this.x0 = fullsys.x0.clone;
                end
                this.C = fullsys.C.clone;
            end
            
            % Check whether projection was setup for this system
            if ~isempty(V)
                this.R = this.compileReconstructionMatrix(V);
                
                % Project input
                if ~isempty(fullsys.B)
                    this.B = fullsys.B.project(V,W);
                end
                % Project linear part
                if ~isempty(fullsys.A)
                    this.A = fullsys.A.project(V,W);
                end
                % Project the approximated CoreFun of the full model if exists
                if ~isempty(fullmodel.Approx)
                    this.f = fullmodel.Approx.project(this.R,W);
                elseif ~isempty(fullsys.f)
                    % Otherwise project the models' full function.
                    this.f = fullsys.f.project(this.R,W);
                end
                % Project mass matrix
                if ~isempty(fullsys.M)
                    this.M = fullsys.M.project(V,W);
                end
                % Project algebraic constraints function
                % This needs only affect the getStateJacobian function of
                % ACoreFun
                if ~isempty(fullsys.g)
                    this.g = fullsys.g.project(this.R,1);
                    this.g.setSystem(this);
                end
            else
                % Only use approximated version if set
                if ~isempty(fullmodel.Approx)
                    this.f = fullmodel.Approx;
                end
            end
            % Set the reduced system as base system for the projected
            % ACoreFun - it will need to change references to the full
            % system itself if required
            this.f.setSystem(this);
            
            % Set the plot-wrapper (uses the plot method from the full
            % system)
            this.plotPtr = @fullsys.plot;
        end
        
        function updateSparsityPattern(this)
            % Reduced systems are not sparse anymore
            this.SparsityPattern = [];
        end
        
        function plot(this, model, t, y)
            % Unless overridden for specific reduced system plots this
            % method just calls the plot method of the original full system
            this.plotPtr(model, t, y);
        end
        
        function dx = ODEFun(this,t,x)
            est = this.Model.ErrorEstimator;
            haveest = ~isempty(est) && est.Enabled;
            if haveest
                xall = x;
                x = x(1:end-est.ExtraODEDims,:);    
            end
            
            dx = ODEFun@models.BaseFirstOrderSystem(this,t,x);
            
            if haveest
                dx = [dx; est.evalODEPart(xall, t, this.u(t))];
            end
        end
        
        function J = getJacobian(this, t, xc)
            est = this.Model.ErrorEstimator;
            haveest = ~isempty(est) && est.Enabled;
            if haveest
                xc = xc(1:end-est.ExtraODEDims,:);    
            end
            J = getJacobian@models.BaseFirstOrderSystem(this, t, xc);
            if haveest
                J = blkdiag(J,0);
            end
        end
        
        function x0 = getX0(this, mu)
            % Gets the initial value of `x_0(\mu)`.
            %
            % If the estimator is enabled, x0 is extended by the e0
            % components of the error estimator.
            
            x0 = getX0@models.BaseFirstOrderSystem(this, mu);
            m = this.Model;
            if ~isempty(m.ErrorEstimator) && m.ErrorEstimator.Enabled
                x0 = [x0; m.ErrorEstimator.getE0(mu)];
            end
        end
        
        function fsys = get.FullSystem(this)
            fsys = this.Model.FullModel.System;
        end
    end
    
    methods(Access=protected)
        
        function R = compileReconstructionMatrix(this, V)
            R = V;
        end
        
        function x0proj = projectx0(~, x0, SV, SW)
            % Overridden in ReducedSecondOrderSystem
            x0proj = x0.project(SV,SW);
        end
        
        function updateDimensions(this)
            rm = this.Model;
            if ~isempty(rm.V)
                % This is the projection case
                this.NumStateDofs = size(rm.V,2);
            else
                % No projection case
                this.NumStateDofs = rm.FullModel.System.NumStateDofs;
            end
            this.NumAlgebraicDofs = rm.FullModel.System.NumAlgebraicDofs;
            this.NumTotalDofs = this.NumStateDofs + this.NumAlgebraicDofs;
        end
        
        function validateModel(this, model)%#ok
            % Overrides the validateModel function in BaseFirstOrderSystem.
            % 
            % No need to call the superclass here as ReducedModel is a
            % child of BaseModel.
            if ~isa(model, 'models.ReducedModel')
                error('The Model property has to be a child of models.ReducedModel');
            end
        end
    end
end

