classdef ReducedSystem < models.BaseDynSystem
    %ReducedSystem: A KerMor reduced dynamical system.
    % 
    % @author Daniel Wirtz @date 17.03.2010
    %
    % @change{0,4,dw,2011-05-13} Removed the old 'getODEFun' function and replaced it by a direct
    % implementation 'ODEFun'. This enables to avoid nested function handles with in turn allow for
    % a speedup of reduced simulations by almost a factor of 2.
    %
    % See also: models.BaseDynSystem
    
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
            
            % Clones the full system's basic (all but functions)
            % properties
            this.Params = fullsys.Params;
            this.Inputs = fullsys.Inputs;
            this.MaxTimestep = fullsys.MaxTimestep;
            this.StateScaling = fullsys.StateScaling;
            
            % Copy component handles (it's NOT cloning them!)
            % This is the default as if no reduction methods are
            % applied (sensless but yet allowed) an identical system
            % will be the result.
            this.f = fullsys.f;
            this.B = fullsys.B;
            this.C = fullsys.C;
            
            % Forwards the x0 evaluation to the original model's x0 function.
            if ~isempty(this.Model.W)
                this.x0 = @(mu)this.Model.W'*fullsys.x0(mu);
            else
                this.x0 = fullsys.x0;
            end
            
            % Set the plot-wrapper (uses the plot method from the full
            % system)
            this.plotPtr = @fullsys.plot;
            
            % Create local workspace copy
            V = this.Model.V;
            W = this.Model.W;
            
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
                this.C = fullsys.C.project(V,W);
                
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
                y = [y; est.evalODEPart(x,t,this.mu,ut)];
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

