classdef ReducedSystem < models.BaseDynSystem
    %ReducedSystem: A KerMor reduced dynamical system.
    %
    % Any state space scaling is included into the reduced model upon creation, as scaling is
    % possibly of the large systems' dimension. See @ref scaling for more information.
    % 
    % @author Daniel Wirtz @date 17.03.2010
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
            
            % The state scaling for the reduced system is one, as the scaling matrices have been
            % incorporated into the V,W matrices via V := SV and W = S^-1W inside the ReducedModel
            % setModel method.
            this.StateScaling = 1;
            s = fullmodel.System.StateScaling;
            d = size(fullmodel.Data.V,1);
            if isscalar(s)
                s(1:d,1) = s;
            end
            % Incorporate the state scaling into the projection matrices.
            % They are used to project the x0 initial values and output C.
            SV = spdiags(s,0,d,d) * V;
            SW = spdiags(1./s,0,d,d) * W;
            
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
            this.C = fullsys.C;
            
            % Forwards the x0 evaluation to the original model's x0 function.
            if ~isempty(W)
                this.x0 = fullsys.x0.project(SV,SW);
            else
                this.x0 = fullsys.x0.clone;
            end
            
            % Set the plot-wrapper (uses the plot method from the full
            % system)
            this.plotPtr = @fullsys.plot;
            
            % Figure whether preojection was setup for this system
            if ~isempty(fullmodel.SpaceReducer)
                if isempty(V) || isempty(W)
                    error(['Model has a SpaceReducer but no projection'...
                        'data is given. Forgot to call offlineGenerations?']);
                end
                % Project input/output
                if ~isempty(fullsys.B)
                    this.B = fullsys.B.project(V,W);
                end
                % We have a C in each case, mostly StdOutputConv.
                this.C = fullsys.C.project(SV,SW);
                
                % Project the approximated CoreFun of the full model if exists
                if ~isempty(fullmodel.Approx)
                    this.f = fullmodel.Approx.project(V,W);
                else
                    % Otherwise at least try to project the models' full
                    % function.
                    this.f = fullsys.f.project(V,W);
                end
            else
                % Only use approximated version if set
                if ~isempty(fullmodel.Approx)
                    this.f = fullmodel.Approx;
                end
            end
        end
                
        function y = ODEFun(this, t, x)
            % Overrides the default implementation in BaseDynSystem and extends the functionality of
            % the ODE function if any error estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est) && est.Enabled
%                 V = 1;
%                 if ~isempty(this.Model.V)
%                     V = this.Model.V;
%                 end
                % Eval f part
                y = this.f.evaluate(x(1:end-est.ExtraODEDims,:),t,this.mu);
                % See if Bu is used
                if ~isempty(this.u)
                    ut = this.u(t);
                    y = y + this.B.evaluate(t,this.mu)*ut;
                else
                    ut = [];
                end
                % Extend by error estimator part
                y = [y; est.evalODEPart(x, t, this.mu, ut)];
            else
                % If no estimator is used or is disabled just call the "normal" ODE function from the
                % base class.
                y = ODEFun@models.BaseDynSystem(this, t, x);
            end
        end
        
        function plot(this, model, t, y)
            % Unless overridden for specific reduced system plots this
            % method just calls the plot method of the original full system
            this.plotPtr(model, t, y);
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

