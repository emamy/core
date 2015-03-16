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
    % @todo maybe make W dependent and return V if W has not been explicitly set
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        % @type matrix<double> @default []
        V;
        
        % The biorthogonal matrix for V, i.e. `W^tV = I_d`
        %
        % If empty, and `V` is set, `W=V` is assumed.
        %
        % @type matrix<double> @default []
        W;
        
        % The originally used parameter samples.
        %
        % The parameter samples that have been used for computation of the
        % reduced model.
        %
        % @type matrix<double>
        ParamSamples;
    end
    
    properties(Dependent)
        % The error estimator for the reduced model
        %
        % @type error.BaseEstimator @default []
        ErrorEstimator;
    end
    
    properties(Access=private)
        % Internal variable for the dependent property ErrorEstimator
        %
        % @type error.BaseEstimator
        fErrorEstimator;
    end
    
    methods
        
        function this = ReducedModel(fullmodel, target_dim)
            this = this@models.BaseModel;
            if nargin > 0
                this.setFullModel(fullmodel, target_dim);
            end
        end
        
        function delete(this)
            this.ErrorEstimator = [];
            this.V = [];
            this.W = [];
            % Delete the full model with it's ModelData instance at last (V, W might have data
            % there)
            this.FullModel = [];
        end
        
    end
    
    methods(Sealed)
        
        function setFullModel(this, fullmodel, target_dim)
            % Creates a new reduced model instance from a full model.
            %
            % Ensure that offlineGenerations have been called on the full
            % model to provide any necessary reduction data.
            %
            % Parameters:
            % fullmodel: A full model where the reduced model shall be
            % created from. @type models.BaseFullModel
            % target_dim: The target dimension `d` of the reduced model.
            % Uses the first `d` columns of the projection matrices `V`
            % (and `W` if set, respectively) to create a subspace-projected
            % reduced model.
            %
            % See also: models.BaseFullModel
            % models.BaseFullModel.buildReducedModel
            % models.BaseFullModel.offlineGenerations            
            if nargin == 0 || ~isa(fullmodel,'models.BaseFullModel')
                error('ReducedModel instances require a full model to construct from.');
            elseif isempty(fullmodel.Approx) && isempty(fullmodel.SpaceReducer)
                error('No reduction methods found on full model. No use in building a reduced model from it.');
            end
            
            fd = fullmodel.Data;
            ns = length(fd.ProjectionSpaces);
            
            %
%             algdofs = fullmodel.System.AlgebraicConditionDoF;
%             numalgdofs = length(algdofs);
            
            disp('Building reduced model...');
            
            
            % If no target dimension(s) are given, assume full size (for
            % each subspace)
            if nargin < 3 || isempty(target_dim)
                for k = 1:ns
                    target_dim(k) = fd.ProjectionSpaces(k).Size;
                end
                % If only a scalar is given, assume a total desired size of
                % target_dim and equally split sizes for each subspace
            elseif isscalar(target_dim)
                target_dim = ones(ns,1)*floor(target_dim/ns);
            end
            for k = 1:ns
                if target_dim(k) > fd.ProjectionSpaces(k).Size
                    warning('Target size %d for subspace %d too large! Using max value of %d',...
                        target_dim(k),k,fd.ProjectionSpaces(k).Size);
                    target_dim(k) = fd.ProjectionSpaces(k).Size;
                end
            end
            
            % IMPORTANT: Assign any model properties that are used during
            % the creation of the reduced system! (i.e. V,W are used in the
            % constructor of the reduced System)
            this.FullModel = fullmodel;
            % Update name ;-)
            this.Name = ['Reduced: ' fullmodel.Name];
            
            %% Copy common values from the full model
            this.T = fullmodel.T;
            this.dt = fullmodel.dt;
            this.tau = fullmodel.tau;
            this.isStatic = fullmodel.isStatic;
            this.DefaultMu = fullmodel.DefaultMu;
            this.DefaultInput = fullmodel.DefaultInput;
            fms = fullmodel.ODESolver;
            if isa(fms,'solvers.SemiImplicitEuler')
                this.ODESolver = solvers.SemiImplicitEuler(this);
            elseif isa(fms,'solvers.FullyImplEuler')
                this.ODESolver = solvers.FullyImplEuler(this);
            else
                this.ODESolver = fms;
            end
            this.G = fullmodel.G;
            
            %% Compute the effective projection matrix
            % As we might have multiple projection subspaces, we need to
            % assemble a global projection matrix (block-like if separate
            % subspace dimensions are "sorted" in full model).
            % Any not-projected DoFs are forwarded as-is by identity. This
            % automatically deals with DAE constraints as they are not set
            % to be reduced anyways.
            V = []; W = []; %#ok<*PROP>
            if ns > 0
                dim = fullmodel.System.StateSpaceDimension;
                V = zeros(dim,sum(target_dim));
                W = V;
                offset = 0;
                done = false(dim,1);
                for k=1:ns
                    s = fd.ProjectionSpaces(k);
                    sel = 1:target_dim(k);
                    V(s.Dimensions,offset + sel) = s.V(:,sel);
                    %if ~isempty(s.W)
                        W(s.Dimensions,offset + sel) = s.W(:,sel);
                    %end
                    offset = offset + target_dim(k);
                    done(s.Dimensions) = true;
                end
                % Insert identity for any remaining dimensions
                unreduced = find(~done);
                nunred = length(unreduced);
                V(unreduced,end+(1:nunred)) = eye(nunred);
                W(unreduced,end+(1:nunred)) = eye(nunred);
            end
            this.V = V;
            this.W = W;
            %this.V = data.FileMatrix(V,'Dir',fd.DataDirectory);
            %this.W = data.FileMatrix(W,'Dir',fd.DataDirectory);
            
            this.ParamSamples = fullmodel.Data.ParamSamples;
            
            %% Create a new reducedSystem passing this reduced model
            this.System = fullmodel.System.buildReducedSystem(this);
            
            % Use the error estimator that has been precomputed in 
            % the full model
            if ~isempty(fullmodel.ErrorEstimator)
                this.ErrorEstimator = fullmodel.ErrorEstimator.prepareForReducedModel(this);
            end
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
            % @type rowvec<double>
            % x: The state variables at the corresponding times t. @type matrix<double>
            % ctime: The time needed for computation. @type double
            
            % Call constat pre-computations
            cpre = 0;
            if ~isempty(this.ErrorEstimator)
                this.ErrorEstimator.clear;
                if this.ErrorEstimator.Enabled
                    cpre = this.ErrorEstimator.prepareConstants(mu, inputidx);
                end
            end
            
            % Call inherited method (actual work)
            [t, xext, ctime] = computeTrajectory@models.BaseModel(this, mu, inputidx);
            
            % Split up results; the last rows of the ode solution contain
            % any online-computable errors
            cpost = 0;
            if ~isempty(this.ErrorEstimator) && this.ErrorEstimator.Enabled
                x = xext(1:end-this.ErrorEstimator.ExtraODEDims,:);
                cpost = this.ErrorEstimator.postProcess(xext, t, inputidx);
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
            
            if ~isempty(this.ErrorEstimator) && this.ErrorEstimator.Enabled
                x0 = [x0; this.ErrorEstimator.getE0(mu)];
            end
        end
        
    end
    
    methods(Sealed)
        function [f, ax] = plot(this, t, y, varargin)
            [f,ax] = this.FullModel.plot(t, y, varargin{:});
        end
        
        function [f, ax] = plotState(this, t, y, varargin)
            [f,ax] = this.FullModel.plotState(t, y, varargin{:});
        end
        
        function [f, ax] = plotSingle(this, t, y, varargin)
            [f,ax] = this.FullModel.plotSingle(t, y, varargin{:});
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
            Utils.saveFigure(h, fullfile(folder,file),'png');
            file = [file '.png'];
            close(h);
        end
        
        function e = get.ErrorEstimator(this)
            e = this.fErrorEstimator;
        end
        
        function set.ErrorEstimator(this, value)
            %
            % Parameters:
            % value: @type error.BaseEstimator
            if ~isempty(value)
                if ~isa(value,'error.BaseEstimator')
                    error('The ErrorEstimator property must be a subclass of the error.BaseEstimator class.');
                end
                msg = value.validModelForEstimator(this.FullModel);
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

