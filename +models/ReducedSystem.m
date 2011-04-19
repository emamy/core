classdef ReducedSystem < models.BaseDynSystem
    %REDUCEDSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
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
            this.x0 = fullsys.x0;
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
                
                % Project the initial value function
                % Dont use this.Data.V but reduced.V since the function
                % handle stores the local workspace and after loading the
                % reduced model the Data field will be [].
                % See ReducedModel.save for details.
                this.x0 = @(mu)W'*fullsys.x0(mu);
            else
                % Only use approximated version if set
                if ~isempty(fullmodel.Approx)
                    this.f = fullmodel.Approx;
                end
            end
        end
        
        function odefun = getODEFun(this, mu, inputidx)
            % Overrides the default implementation in BaseDynSystem and
            % extends the functionality of the ODE function if any error
            % estimators are enabled.
            est = this.Model.ErrorEstimator;
            if ~isempty(est)
                if est.Enabled
                    xdims = est.ExtraODEDims;
                    % System without inputs
                    if this.InputCount == 0 || isempty(this.B)
                        odefun = @(t,x)([this.f.evaluate(x(1:end-xdims,:),t,mu); est.evalODEPart(x,t,mu)]);
                    else
                        % generates the ode function for given parameter and input function
                        u = this.Inputs{inputidx};
                        odefun = @(t,x)([this.f.evaluate(x(1:end-xdims,:),t,mu) + this.B.evaluate(t,mu)*u(t); est.evalODEPart(x,t,mu,u(t))]);
                    end
                    return;
                else
                    % Clear the last error
                    est.clear;
                end
            end
            % If no estimator is used or is disabled just call the "normal" ODE
            % function generator from the base class.
            odefun = getODEFun@models.BaseDynSystem(this, mu, inputidx);
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
    
%     methods(Static,Access=protected)
%         function obj = loadobj(s)
%             obj = models.ReducedSystem;
%             ALoadable.getObjectDict(s.ID) = obj;
%             % Load BaseModel's properties
%             obj = models.BaseDynSystem.loadobj(s, obj);
%             % Load local properties
%             obj.plotPtr = s.plotPtr;
%         end
%     end
end

