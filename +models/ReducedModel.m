classdef ReducedModel < models.BaseModel
    % The KerMor reduced model class
    %
    % In difference to the BaseModels this class cannot be extended. All
    % needed functionality for reduced models is included here (simulation,
    % error computation etc) and no specific subclasses for other models
    % need to be implemented. One can obtain a reduced model from a full
    % model by calling the full model's @code buildReducedModel @endcode
    % function. Dont forget to call @code offlineGenerations @endcode at
    % least once before reduction.
    %
    % See also: BaseModel BaseFullModel
    %
    % @author Daniel Wirtz @date 23.03.2010
    
    properties(SetAccess=private)
        
        % The full model this reduced model was created from.
        %
        % Once an instance of a reduced model gets saved to disk, the
        % FullModels .Data and .Approx properties get set to [] for disk
        % space reduction since these properties are not necessarily used
        % by the reduced system. So far, only error computations will take
        % longer since the initial snapshots are not available anymore.
        FullModel;
        
        % The matrix that has been used for projection
        %
        V;
        
        % The biorthogonal matrix for V, i.e. `W^tV = I_d`
        %
        W;
        
        % The originally used parameter samples.
        %
        % The parameter samples that have been used for computation of the
        % reduced model.
        ParamSamples;
    end
    
    properties
        % The error estimator for the reduced model
        %
        %
        %
        % @type error.BaseEstimator
        ErrorEstimator;
    end
    
    methods(Sealed)
        
        function this = ReducedModel(fullmodel)
            if nargin == 1
                this.setFullModel(fullmodel);
            end
        end
        
        function setFullModel(this, fullmodel)
            % Creates a reduced model from a given full model.
            %
            % @docupdate
            if nargin == 0 || ~isa(fullmodel,'models.BaseFullModel')
                error('ReducedModel instances require a full model to construct from.');
            end
            
            % Check if a reduction is available at all
            if isempty(fullmodel.Approx) && isempty(fullmodel.SpaceReducer)
                warning('KerMor:Reducing:noEffectiveReduction',...
                    ['Model setup lacks reduction methods.\n'...
                    'Reduced model will equal the full model (computationally).']);
            end
            
            disp('Start building reduced model...');
            % IMPORTANT: Assign any model properties that are used during
            % the creation of the reduced system! (i.e. V,W are used in the
            % constructor of the reduced System)
            this.FullModel = fullmodel;
            % Update name ;-)
            this.Name = ['Reduced: ' fullmodel.Name];
            
            % Copy common values from the full model
            this.T = fullmodel.T;
            this.dt = fullmodel.dt;
            this.Verbose = fullmodel.Verbose;
            this.ODESolver = fullmodel.ODESolver;
            this.G = fullmodel.G;
            
            % Copy data that is also needed in the reduced model
            this.V = fullmodel.Data.V;
            this.W = fullmodel.Data.W;
            this.ParamSamples = fullmodel.Data.ParamSamples;
            
            % Create a new reducedSystem passing the full and so far
            % initialized reduced models
            this.System = models.ReducedSystem(fullmodel,this);
            
            % Obtain an error estimator. The static method looks for a
            % suitable error estimator or returns the default (=expensive!)
            % estimator.
            this.ErrorEstimator = error.BaseEstimator.getEstimator(this);
        end
        
        function [t,x] = computeTrajectory(this, mu, inputidx)
            % Call parent method for actual work
            % @docupdate
            
            % Clear possibly old data in error estimators
            this.ErrorEstimator.clear;
            
            % Call inherited method (actual work)
            [t,xtmp] = computeTrajectory@models.BaseModel(this, mu, inputidx);
            
            % Split up results; the last row of the ode solution contains
            % any online-computable errors
            if this.ErrorEstimator.Enabled
                x = xtmp(1:end-this.ErrorEstimator.ExtraODEDims,:);
                this.ErrorEstimator.process(t, xtmp, mu, inputidx);
            else
                x = xtmp;
            end
        end
        
        function exo = getExo(this, mu)
            % Computes the norm of the initial error 
            % `E_{y_0}(\mu) = ||C(0,\mu)(I-VW^t)x_0(\mu)||`
            if ~isempty(this.V) && ~isempty(this.W)
                x0 = this.FullModel.System.x0(mu);
                Vy = this.V*(this.W'*x0);
                exo = x0'*this.G*x0 - 2*x0'*this.G*Vy + Vy'*this.G'*Vy;
                exo = sqrt(max(exo,0));
%                 exo = sqrt(sum((C*(x0-this.V*this.W'*x0)).^2));
            else
                exo = 0;
            end
        end
        
%         function [t,e,est] = getError(this, varargin)
%             % Gets the reduced model's error for a specified parameter and
%             % input.
%             %
%             % @todo possibly move into BaseEstimator?
%             
%             [t,x,xr] = this.getTrajectories(varargin{:});
%             % Compute L^2-norm of difference for each timestep
%             e = sqrt(sum((x - xr).^2,1));
%             est = this.ErrorEstimator.LastError;
%         end
        
%         function [t,erel,estrel] = getRelativeError(this, varargin)
%             [t,x,xr] = this.getTrajectories(varargin{:});
%             e = sqrt(sum((x - xr).^2,1));
%             est = this.ErrorEstimator.LastError;
%             x = sqrt(sum(x.^2,1));
%             erel = e./x;
%             estrel = est./x;
%         end        
        
    end
    
    methods(Access=protected,Sealed)
        
        function x0 = getX0(this, mu)
            % Gets the initial value of `x_0(\mu)`.
            %
            % If the estimator is enabled, x0 is extended by the e0
            % components of the error estimator.
            x0 = getX0@models.BaseModel(this, mu);
            if this.ErrorEstimator.Enabled
                x0 = [x0; this.ErrorEstimator.getE0(mu)];
            end
        end
        
    end
    
    methods
        
        function [t,x,xr,time,timer,t_noerr] = getTrajectories(this, mu, inputidx)
            % Debug Method. Computes the trajectories of the full, reduced
            % and reduced without error estimation systems.
            if nargin < 3
                inputidx=[];
                if nargin < 2
                    mu = [];
                end
            end
            % Perform full simulation (computed trajectories are cached!)
            tic;
            [t,x] = this.FullModel.computeTrajectory(mu, inputidx);
            time = toc;
            % Perform reduced simulation
            this.ErrorEstimator.Enabled = false;
            tic;
            this.computeTrajectory(mu, inputidx);
            %xr = this.V*z;%#ok
            t_noerr = toc;
            this.ErrorEstimator.Enabled = true;
            tic;
            [t,z] = this.computeTrajectory(mu, inputidx);
            xr = this.V*z;
            timer = toc;
        end
        
        function set.ErrorEstimator(this, value)
            if ~isa(value,'error.BaseEstimator')
                error('The ErrorEstimator property must be a subclass of the error.BaseEstimator class.');
            end
            msg = value.validModelForEstimator(this);
            if ~isempty(msg)
                error(msg);
            end
            this.ErrorEstimator = value;
        end
    end
    
    % Save & Load
    methods
        function s = saveobj(this)
            % Saves the reduced model to a struct.
            % For the save process of the reduced model the full model's
            % Data (=ModelData) and Approx properties are not needed. This
            % is the fastest way to ensure that the reduced model can still
            % have access to all important features of the full model but
            % uses less disk space.
            
            s.FullModel = this.FullModel;
            s.FullModel = this.FullModel;
            s.V = this.V; r.W = this.W;
            s.ParamSamples = this.ParamSamples;
            s.ErrorEstimator = this.ErrorEstimator;
        end
        
        function this = reload(this,s)
            this.FullModel = s.FullModel;
            this.V = s.V; r.W = s.W;
            this.ParamSamples = s.ParamSamples;
            this.ErrorEstimator = s.ErrorEstimator;
        end
    end
    
    methods (Static)
        function obj = loadobj(S)
            obj = reload(models.ReducedModel,S);
        end
    end
    
end
