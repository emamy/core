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
    %
    % @new{0,5,dw,2011-08-23} Added a createImage method that per default
    % saves the KerMor Logo with the models name in the
    % KerMor.TempDirectory. This is used for example in general.AppExport
    % to create an representing image for a reduced model.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing     
    
    properties(SetAccess=private)
        
        % The full model this reduced model was created from.
        %
        % Once an instance of a reduced model gets saved to disk, the
        % FullModels .Data and .Approx properties get set to [] for disk
        % space reduction since these properties are not necessarily used
        % by the reduced system. So far, only error computations will take
        % longer since the initial snapshots are not available anymore.
        %
        % @type models.BaseFullModel
        FullModel;
        
        % The matrix that has been used for projection
        %
        % @type matrix
        V;
        
        % The biorthogonal matrix for V, i.e. `W^tV = I_d`
        %
        % @type matrix
        W;
        
        % The originally used parameter samples.
        %
        % The parameter samples that have been used for computation of the
        % reduced model.
        %
        % @type matrix
        ParamSamples;
    end
    
    properties(Dependent)
        % The error estimator for the reduced model
        %
        % @type error.BaseEstimator
        ErrorEstimator;
    end
    
    properties(Access=private)
        fErrorEstimator;
    end
    
    methods(Sealed)
        
        function this = ReducedModel(fullmodel)
            % Creates a new reduced model instance.
            %
            % Optionally, a models.BaseFullModel subclass can be passed to
            % create this new reduced model from.
            %
            % Parameters:
            % fullmodel: A full model where the reduced model shall be
            % created from. [Optional]
            this = this@models.BaseModel;
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
            this.tau = fullmodel.tau;
            this.ODESolver = fullmodel.ODESolver;
            this.G = fullmodel.G;
            
            % Copy data that is also needed in the reduced model
            this.V = fullmodel.Data.V;
            this.W = fullmodel.Data.W;
            
            this.ParamSamples = fullmodel.Data.ParamSamples;
            
            % Create a new reducedSystem passing this reduced model
            this.System = models.ReducedSystem(this);
            
            % Obtain an error estimator. The static method looks for a
            % suitable error estimator or returns the default (=expensive!)
            % estimator.
            this.ErrorEstimator = error.BaseEstimator.getEstimator(this);
        end
        
        function [t, x, ctime] = computeTrajectory(this, mu, inputidx)
            % Computes an approximated solution/trajectory for the given mu and inputidx in the
            % coefficient space.
            %
            % Parameters:
            % mu: The parameter `\mu` for the simulation @type colvec
            % inputidx: The integer index of the input function to use. If
            % more than one inputs are specified this is a necessary argument. @type integer
            %
            % Return values:
            % t: The times at which the model was evaluated. Will equal the property Times
            % @type rowvec
            % x: The state variables at the corresponding times t. @type matrix<double>
            % ctime: The time needed for computation. @type double
            
            % Clear possibly old data in error estimators
            this.ErrorEstimator.clear;
            
            % Call constat pre-computations
            cpre = this.ErrorEstimator.prepareConstants(mu, inputidx);
            
            % Call inherited method (actual work)
            [t, xext, ctime] = computeTrajectory@models.BaseModel(this, mu, inputidx);
            
            % Split up results; the last rows of the ode solution contain
            % any online-computable errors
            cpost = 0;
            if this.ErrorEstimator.Enabled
                x = xext(1:end-this.ErrorEstimator.ExtraODEDims,:);
                cpost = this.ErrorEstimator.postProcess(xext, t, mu, inputidx);
            else
                x = xext;
            end
            ctime = ctime + cpre + cpost;
        end
             
        function saveFinal(this, filename)
            % Saves this reduced model for final use.
            %
            % Using this method produces a small-sized file only suitable for online simulations.
            % Most other functionalities or analysis will not work anymore as all large data has
            % been discarded during the saving process.
            a = this.FullModel.Approx;
            d = this.FullModel.Data;
            
            this.FullModel.Approx = [];
            this.FullModel.Data = [];
            
            save(filename, 'this');
            
            this.FullModel.Approx = a;
            this.FullModel.Data = d;
        end
    end
    
    methods(Access=protected,Sealed)
        
        function x0 = getX0(this, mu)
            % Gets the initial value of `x_0(\mu)`.
            %
            % If the estimator is enabled, x0 is extended by the e0
            % components of the error estimator.
            
            x0 = getX0@models.BaseModel(this, mu);
            %x0 = this.System.x0.evaluate(mu);
            if this.ErrorEstimator.Enabled
                x0 = [x0; this.ErrorEstimator.getE0(mu)];
            end
        end
        
    end
    
    methods(Sealed)
        function plot(this, t, y)
            this.FullModel.plot(t, y);
        end
        
        function plotSingle(this, t, y, varargin)
            this.FullModel.plotSingle(t, y, varargin{:});
        end
    end
    
    methods
        
        function [file, folder] = createImage(this)
            % This method can be invoked to obtain an image for the current
            % model.
            %
            % If not overridden in subclasses, this method creates an png image
            % using the KerMor logo and the model's name as text in it.
            
            n = this.Name;
            h = KerMorLogo;
            text(-6,5,1.05,n,'FontSize',14,'FontWeight','bold');
            folder = KerMor.App.TempDirectory;
            file = 'modelimage';
            general.Utils.saveFigure(h, fullfile(folder,file),'png');
            file = [file '.png'];
            close(h);
        end
        
        function e = get.ErrorEstimator(this)
            e = this.fErrorEstimator;
        end
        
        function set.ErrorEstimator(this, value)
            if ~isempty(value)
                if ~isa(value,'error.BaseEstimator')
                    error('The ErrorEstimator property must be a subclass of the error.BaseEstimator class.');
                end
                msg = value.validModelForEstimator(this);
                if ~isempty(msg)
                    error(msg);
                end
            end
            this.fErrorEstimator = value;
        end
    end
    
%     methods(Static,Access=protected)
%         function obj = loadobj(s)
%             % Creates a new reduced model instance and loads its properties
%             % from the given struct. 
%             %
%             % This method is only implemented as the property assignment
%             % order is important for a reduced model (example: cannot set
%             % an error estimator before having set the models system)
%             obj = models.ReducedModel;
%             
%             % Load BaseModel's properties
%             obj = loadobj@models.BaseModel(s, obj);
%             
%             % Load local properties
%             di = ALoadable.getObjectDict;
%             ex = di(s.FullModel.ID);
%             if ~isempty(ex)
%                 obj.FullModel = ex;
%             else
%                 obj.FullModel = s.FullModel;
%             end
%             obj.V = s.V; obj.W = s.W;
%             obj.ParamSamples = s.ParamSamples;
%             obj.ErrorEstimator = s.ErrorEstimator;
%         end
%     end
    
end

