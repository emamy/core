classdef ReducedSystem < models.BaseDynSystem
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
% See also: models.BaseDynSystem
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing    
    
    properties(Access=private)
        plotPtr;
    end
    
    methods
        function this = ReducedSystem(rmodel)
            % Creates a new ReducedSystem instance.
            %
            % Parameters:
            % rmodel: [Optional] The reduced model to create the reduced
            % system from.
            this = this@models.BaseDynSystem(rmodel);
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
        
        function setReducedModel(this, rmodel)
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
            fullmodel = rmodel.FullModel;
            
            fullsys = fullmodel.System;
            
            % Create local workspace copy (if pointers are used, dont store model.V in it..)
            V = this.Model.V;
            W = this.Model.W;
            
            if ~isempty(fullmodel.SpaceReducer) && (isempty(V) || isempty(W))
                error(['Model has a SpaceReducer but no projection'...
                        'data V,W is given. Forgot to call offlineGenerations?']);
            end
            
            % The state scaling for the reduced system is one, as the scaling matrices have been
            % incorporated into the V,W matrices via V := SV and W = S^-1W inside the ReducedModel
            % setModel method.
            this.StateScaling = 1;
            SV = 1; SW = 1; % Default: no scaling
            s = fullmodel.System.StateScaling;
            if s ~= 1
                if isscalar(s)
                    dim = fullsys.x0.evaluate(fullsys.getRandomParam);
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
            if SV ~= 1 
                this.x0 = fullsys.x0.project(SV,SW);
                this.C = fullsys.C.project(SV,SW);
            else
                this.x0 = fullsys.x0.clone;
                this.C = fullsys.C.clone;
            end
            
            % Check whether projection was setup for this system
            if ~isempty(V)
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
                    this.f = fullmodel.Approx.project(V,W);
                else
                    % Otherwise project the models' full function.
                    this.f = fullsys.f.project(V,W);
                end
                % Project mass matrix
                if ~isempty(fullsys.M)
                    this.M = fullsys.M.project(V,W);
                end
            else
                % Only use approximated version if set
                if ~isempty(fullmodel.Approx)
                    this.f = fullmodel.Approx;
                end
            end
            
            % Set the plot-wrapper (uses the plot method from the full
            % system)
            this.plotPtr = @fullsys.plot;
        end
        
        function plot(this, model, t, y)
            % Unless overridden for specific reduced system plots this
            % method just calls the plot method of the original full system
            this.plotPtr(model, t, y);
        end
                
        function y = ODEFun_A(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
                y = [this.A.evaluate(x(1:end-est.ExtraODEDims,:), t, this.mu);...
                     est.evalODEPart(x, t, this.mu, [])];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun_A@models.BaseDynSystem(this, t, x);
            end
        end
        
        function y = ODEFun_f(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
                y = [this.f.evaluate(x(1:end-est.ExtraODEDims,:),t,this.mu);...
                     est.evalODEPart(x, t, this.mu, [])];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun_f@models.BaseDynSystem(this, t, x);
            end
        end
        
        function y = ODEFun_Af(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
                y = [this.A.evaluate(x(1:end-est.ExtraODEDims,:), t, this.mu) + ...
                     this.f.evaluate(x(1:end-est.ExtraODEDims,:), t, this.mu); ...
                     est.evalODEPart(x, t, this.mu, [])];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun_Af@models.BaseDynSystem(this, t, x);
            end
        end
        
        function y = ODEFun_AB(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
                ut = this.u(t);    
                y = [this.A.evaluate(x(1:end-est.ExtraODEDims,:), t, this.mu) + ...
                     this.B.evaluate(t,this.mu)*ut; ...
                     est.evalODEPart(x, t, this.mu, ut)];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun_AB@models.BaseDynSystem(this, t, x);
            end
        end
        
        function y = ODEFun_fB(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
                ut = this.u(t);    
                y = [this.f.evaluate(x(1:end-est.ExtraODEDims,:),t,this.mu) + ...
                     this.B.evaluate(t,this.mu)*ut; ...
                     est.evalODEPart(x, t, this.mu, ut)];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun_fB@models.BaseDynSystem(this, t, x);
            end
        end
        
        function y = ODEFun_AfB(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
                ut = this.u(t);    
                y = [this.A.evaluate(x(1:end-est.ExtraODEDims,:), t, this.mu) + ...
                     this.f.evaluate(x(1:end-est.ExtraODEDims,:), t, this.mu) + ...
                     this.B.evaluate(t,this.mu)*ut; ...
                     est.evalODEPart(x, t, this.mu, ut)];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun_AfB@models.BaseDynSystem(this, t, x);
            end
        end
    end
    
    methods(Access=protected, Sealed)
        function validateModel(this, model)%#ok
            % Overrides the validateModel function in BaseDynSystem.
            % 
            % No need to call the superclass here as ReducedModel is a
            % child of BaseModel.
            if ~isa(model, 'models.ReducedModel')
                error('The Model property has to be a child of models.ReducedModel');
            end
        end
    end
end

