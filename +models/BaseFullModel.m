classdef BaseFullModel < models.BaseModel & IParallelizable
% The base class for any KerMor detailed model
%
% Implementers of custom models are to inherit from this base class
% in order for it to work with KerMor.
% For custom models, the properties of this class (combined with
% those from BaseModel) can be set to influence the model behaviour
% and reduction methods.
% For the implementation of custom dynamical systems, refer to
% BaseDynSystem.
%
% Setting default values for properties that are handle classes will have the negative
% side-effect of having each instance of BaseFullModel initialized with the SAME instance of the
% Sampler, Approx etc. instances which of course is NOT desireable.
% See http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_oop/bsdtcz7.html#bsdu1g9-1
% for details.
%
% @todo build in time-tracking for offline phase etc
%
% @todo If enabletrajectorycaching is off, the model complains about no
% finding file system folders for trajectories etc after re-loading
%
% See also: models BaseModel BaseDynSystem
%
% @author Daniel Wirtz @date 16.03.2010
%
% @new{0,7,dw,2013-05-28} Added an AutoSave flag for intermediate saving during offline phases.
%
% @new{0,6,dw,2012-10-10} Added the SaveTag property on this level. It is automatically set
% according to the model name whenever it is still empty and the models name is set.
%
% @new{0,6,dw,2012-10-08} Added a new read-only property models.BaseFullModel.OfflinePhaseTimes
% that contains the last offline computation times in a row vector.
%
% @change{0,6,dw,2012-04-26} Using ProcessIndicator for training data
% generation now.
%
% @change{0,6,dw,2011-12-14} Now listening for changes to T and dt of the model and clearing
% the model training data if either value changes
%
% @new{0,5,dw,2011-10-15} Added a new method getTrajApproxError that computes the approximation
% error of the Approx class against the full trajectory. This method provides generic means to
% assess the approximation quality of approx.BaseApprox classes for full trajectories.
%
% @change{0,5,dw,2011-10-16} Fixed the parallel computation of
% BaseFullModel.off2_genTrainingData so that it also works with
% data.FileTrajectoryData (parallel execution did not sync the hashmaps, now
% running data.FileTrajectoryData.consolidate fixes this)
%
% @change{0,5,dw,2011-08-04} Adopted the off2_genTrainingData method to the new data.ATrajectoryData
% structure. Now all trajectories are stored either in memory or disk, and the data.ATrajectoryData
% classes take care of storage.
%
% @new{0,4,dw,2011-05-31} Added the models.BaseFullModel.OfflinePhaseTimes property.
%
% @change{0,4,sa,2011-05-11} Implemented setters for the
% preApproximationTrainingCallback and
% postApproximationTrainingCallback
%
% @new{0,4,dw,2011-05-06} Small improvements to the DPCS, the correct links to the properties
% defining classes are now used. Further a link to the correct class is created if the property
% wsa inherited from a superclass.
%
% @new{0,4,dw,2011-05-04} Added a new property models.BaseFullModel.TrainingInputs. Now one can
% specifiy which defined inputs are to be used within offline generations.
%
% @change{0,3,dw,2011-04-26} The property changed descriptions now contain links to the
% respective classes containing the property.
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @new{0,3,dw,2011-04-21} - Implemented the recursive property changed search on a
% per-BaseFullModel basis as made possible by the changes in KerMorObject.
% - Added a new method models.BaseFullModel.printPropertyChangedReport that gives detailed
% information about any properties of the model and its subclasses that remained set to their
% default value since initialization.
% - Overloaded the models.BaseModel.computeTrajectory method in order to previously run the
% property changed checks
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing    
%
% @todo put addlistener methods for T,dt change into loadobj!
%
% @todo implement callbacks for the ODE solvers that automatically set computed f-values on
% simulation points if desired
    
    properties
        % The full model's data container.
        % Defaults to an empty container.
        %
        % @propclass{data} The model's data container
        %
        % @default data.ModelData @type data.ModelData
        %
        % See also: data.ModelData
        Data;
        
        % Export instance for possible export of this model to JKerMor.
        %
        % Leave empty if no export is possible.
        %
        % @propclass{data} Only used if model can be exported to JaRMoS
        %
        % @type export.JKerMorExport @default []
        JKerMorExport;
    end
    
    properties(SetObservable)
        % The sampling strategy the Model uses
        %
        % @propclass{important} The sampling strategy affects the quality
        % of the parameter samples used for training.
        %
        % @default sampling.RandomSampler @type sampling.BaseSampler
        %
        % See also: sampling BaseSampler
        Sampler;
        
        % The reduction algorithm for subspaces
        %
        % @propclass{critical} The subspace computation strategy is
        % critical for the quality of the projection subspace
        %
        % @default spacereduction.PODGreedy
        % @type spacereduction.BaseSpaceReducer
        %
        % See also: spacereduction spacereduction.BaseSpaceReducer
        SpaceReducer;
        
        % The approximation method for the CoreFunction
        %
        % @propclass{critical} The approximation technique is critical for
        % the quality of the reduced model.
        %
        % @default [] @type approx.BaseApprox
        %
        % See also: approx approx.algorithms data.selection
        Approx = [];
        
        % Advanced property. 
        % Must be a function handle taking the current model instance.
        %
        % The handle is invoked before the ApproxfValues and approximation
        % data is computed within off4_computeApproximation.
        %
        % This property allows to manually influence the approximation
        % training data, for example if the training data has repeating
        % entries if the problem at hand is homogeneous over a specific
        % dimension. An example is the models.rbmatlab.RiemannBurgers
        % model, which is in fact homogeneous in y-direction. Thus, only
        % distinct x-direction points have to be trained for.
        %
        % @propclass{optional} A callback that can be set to execute user
        % defined code before training the approximation.
        %
        % @default []
        % @type function_handle
        %
        % @todo (optional) create class events from this
        preApproximationTrainingCallback;
        
        % Advanced property. 
        % Must be a function handle taking the current model instance.
        %
        % The handle is invoked at the end of the off4_computeApproximation
        % method.
        %
        % @propclass{optional} A callback that can be set to execute user
        % defined code after training the approximation.
        %
        % @default []
        % @type function_handle
        %
        % See also: preApproximationTrainingCallback
        postApproximationTrainingCallback;
        
        % Flag that enables caching of computed trajectories in a simulation cache stored in
        % KerMor's TempDirectory folder.
        %
        % This flag can speedup re-done simulations significantly. However,
        % thus far the cache only looks up the total simulation time,
        % timestep, parameter and input index (if given). Thus, not in all
        % situations a correct recognition of changes is guaranteed (e.g.
        % changing the input function handle will cause the same old
        % trajectory to be loaded).
        % Consequently, this setting is ''disabled'' by default.
        %
        % @propclass{optional} Caching can speedup computations but may
        % require extra disk space.
        %
        % @type logical @default false
        EnableTrajectoryCaching = false;
        
        % The associated error estimator for this model.
        %
        % It is computed automatically, using the order/selection as in
        % error.BaseEstimator.getEstimator
        %
        % @propclass{optional} Using error estimators with reduced models
        % can facilitate reduction quality measurement and can also be
        % needed by some adaptive algorithms.
        %
        % @type error.BaseEstimator @default []
        ErrorEstimator = [];
        
        % Flag that determines whether fxi data should be computed along with the trajectories
        % or not
        %
        % @todo write setter
        %
        % @type logical @default false
        ComputeTrajectoryFxiData = false;
        
        % Flag to enable automatic saving of the model after each individual offline phase
        % step and at other locations prone to data loss due to lengthy computations.
        %
        % Note that once a model has been saved, any file system folders created with that
        % model will persist until manually deleted.
        %
        % @propclass{optional} This setting just decreases the chance of data loss due to
        % intermediate saves.
        %
        % @type logical @default false
        AutoSave = false;
        
        % A custom tag that can be used as a prefix to files for corresponding model
        % identification.
        %
        % If the property models.BaseModel.Name is set and the SaveTag is
        % still empty, a lowercase, trimmed and
        % spaces-replaced-by-underscore SaveTag will be automatically set.
        %
        % @propclass{optional}
        %
        % @type char @default ''
        SaveTag = '';
    end
    
    properties(SetObservable, Dependent)
        % The indices of inputs to use for training data generation.
        %
        % @propclass{optional} Set subindices to restrict offline generation phase to a subset of
        % inputs.
        %
        % @default 1:InputCount
        % @type integer
        TrainingInputs;
    end
    
    properties(SetAccess=private, Dependent)
        % Gets the number of inputs used for training.
        %
        % @type integer
        TrainingInputCount;
    end
    
    properties(SetAccess=private)
        % The computation times for all phases of the last offline generation.
        %
        % Contains a `1\times 6` row vector with the corresponding times for the 'off1_' to 'off_6'
        % phases.
        %
        % @type rowvec<double> @default []
        OfflinePhaseTimes = [];
    end
    
    properties(Access=private)
        fTrainingInputs = [];
        ftcold = struct;
    end
           
    methods
        
        function this = BaseFullModel(name)
            % Creates a new instance of a full model.
            %
            % Parameters:
            % name: The name of the model. @type char
            % @default ''
            this = this@models.BaseModel;
            
            % Register name change listener (might not have the name set
            % already now)
            this.addlistener('Name','PostSet', @this.nameChanged);
            
            % Set name if given
            if nargin == 1
                this.Name = name;
            end
                       
            % Setup default values for the full model's components
            this.Sampler = sampling.RandomSampler;

            this.SpaceReducer = spacereduction.PODGreedy;
            
            this.Data = data.ModelData(this);
            
            % Set default dynamical system
            this.System = models.BaseDynSystem(this);
            
            % Register default properties
            this.registerProps('Sampler','SpaceReducer','Approx',...
                'preApproximationTrainingCallback','postApproximationTrainingCallback',...
                'Data','TrainingInputs','AutoSave');
            
            % Create a listener for the 
            this.addlistener('T','PreSet', @this.intTimeChanging);
            this.addlistener('dt','PreSet', @this.intTimeChanging);
            this.addlistener('T','PostSet', @this.intTimeChanged);
            this.addlistener('dt','PostSet', @this.intTimeChanged);
        end
        
        function delete(this)
            this.Sampler = [];
            this.SpaceReducer = [];
            this.Approx = [];
            this.System = [];
            % Clean up data at last (other instances might still have data inside DataDirectory
            % folder)
            this.Data = [];
            
            delete@models.BaseModel(this);
        end
        
        function off1_createParamSamples(this)
            % Offline phase 1: Sample generation.
            
            % Sampling
            time = tic;
            if ~isempty(this.Sampler)
                this.Data.ParamSamples = this.Sampler.generateSamples(this);
            elseif this.System.ParamCount > 0
                error('No parameter sampler set for a parameterized system. See package sampling and the Sampler property of the BaseFullModel class.');
            end
            this.OfflinePhaseTimes(1) = toc(time);
            this.autoSave;
        end
        
        function off2_genTrainingData(this)
            % Offline phase 2: Snapshot generation for subspace computation
            %
            % @todo
            % - Optimize snapshot array augmentation by preallocation
            % (later will be some storage class)
            % - Remove waitbar and connect to messaging system
            time = tic;
            
            %             if this.Data.SampleCount == 0 && this.TrainingInputCount == 0
            %                 fprintf('BaseFullModel.genTrainingData: No parameters or inputs configured for training, nothing to do.\n');
            %                 time = toc(time);
            %                 return;
            %             end
            
            num_s = max(1,this.Data.SampleCount);
            num_in = max(1,this.TrainingInputCount);
            
            %             if num_s*num_in < matlabpool('size')
            %                 % @\todo: switch back if changed!
            %                 this.ComputeParallel = false;
            %             end
            
            % Clear old trajectory data.
            this.Data.TrajectoryData.clearTrajectories;
            this.RealTimePlotting = false;
            
            %% Parallel - computation
            if this.ComputeParallel
                if this.ComputeTrajectoryFxiData
                    error('ComputeTrajectoryFxiData not yet implemented for parallel execution');
                end
                idxmat = Utils.createCombinations(1:num_s,1:num_in);
                
                fprintf('Starting parallel projection training data computation of %d trajectories on %d workers...\n',size(idxmat,2),matlabpool('size'));
                % Iterate through all param/input combinations
                paridx = idxmat(1,:); % Use sliced variables here
                inidx = idxmat(2,:);
                clear idxmat;
                
                % Notation switch to be reminded that the remote variable will not
                % exist in the local context outside the parfor!
                remote = this;
                n = size(paridx,2);
                parfor idx = 1:n
                    % Check for parameters
                    mu = [];
                    if remote.Data.SampleCount > 0 %#ok<*PFBNS>
                        mu = remote.Data.getParams(paridx(idx));
                    end
                    % Check for inputs
                    inputidx = [];
                    if remote.TrainingInputCount > 0
                        inputidx = inidx(idx);
                    end
                    
                    % Get trajectory
                    [~, x, ctime] = remote.computeTrajectory(mu, inputidx);
                    
                    % Assign snapshot values
                    remote.Data.TrajectoryData.addTrajectory(x, mu, inputidx, ctime);
                end
                % Build the hash map for the local FileTrajectoryData (parfor
                % add them to remotely instantiated FileTrajectoryDatas and does
                % not sync them (not generically possible)
                if isa(this.Data.TrajectoryData,'data.FileTrajectoryData')
                    this.Data.TrajectoryData.consolidate(this);
                end
                
                %% Non-parallel computation
            else
                % Assume no parameters or inputs
                mu = [];
                inputidx = [];
                
                if KerMor.App.Verbose > 0
                    pi = ProcessIndicator(sprintf('Generating projection training data (%d trajectories)...',num_in*num_s),num_in*num_s);
                end
                % Iterate through all input functions
                for inidx = 1:num_in
                    % Iterate through all parameter samples
                    for pidx = 1:num_s 
                        % Check for parameters
                        if this.Data.SampleCount > 0
                            mu = this.Data.getParams(pidx);
                        end
                        % Check for inputs
                        if this.TrainingInputCount > 0
                            inputidx = this.TrainingInputs(inidx);
                        end
                        
                        % Get trajectory (disable caching as trajectories go to this.Data.TrajectoryData
                        % anyways)
                        oldv = this.EnableTrajectoryCaching;
                        this.EnableTrajectoryCaching = false;
                        [t, x, ctime, cache] = this.computeTrajectory(mu, inputidx);
                        this.EnableTrajectoryCaching = oldv;
                        
                        % Store snapshot values, only if not already from trajectory data
                        if cache ~= 1
                            this.Data.TrajectoryData.addTrajectory(x, mu, inputidx, ctime);
                        end
                        
                        % Evaluate fxi on current values if required
                        if this.ComputeTrajectoryFxiData
                            hlp = tic;
                            fx = this.System.f.evaluate(x,t,repmat(mu,1,length(t)));
                            ctime = toc(hlp);
                            this.Data.TrajectoryFxiData.addTrajectory(fx, mu, inputidx, ctime);
                        end
                        
                        if KerMor.App.Verbose > 0
                            pi.step;
                        end
                    end
                end
                if KerMor.App.Verbose > 0
                    pi.stop;
                end
            end
            
            % Input space span computation (if spacereduction is used)
            if ~isempty(this.SpaceReducer) && this.SpaceReducer.IncludeBSpan
                if KerMor.App.Verbose > 0
                    fprintf('Computing input space span...');
                end
                B = this.System.B;
                this.Data.InputSpaceSpan = [];
                if isempty(B)
                    warning('KerMor:BaseFullModel','Cannot compute B span, no System.B set.');
                    return;
                end
                if isa(B,'dscomponents.LinearInputConv')
                    o = general.Orthonormalizer;
                    this.Data.InputSpaceSpan = o.orthonormalize(B.B);
                elseif isa(B,'dscomponents.AffLinInputConv')
                    o = general.Orthonormalizer;
                    tmp = [];
                    for k = 1:B.N
                        tmp = [tmp B.getMatrix(k)];%#ok
                    end
                    this.Data.InputSpaceSpan = o.orthonormalize(tmp);
                else
                    warning('KerMor:BaseFullModel','B span computation: Case %s not yet implemented',class(B));
                end
                if KerMor.App.Verbose > 0
                    fprintf('dimension %d.\n',size(this.Data.InputSpaceSpan,2));
                end
            end
            
            this.OfflinePhaseTimes(2) = toc(time);
            this.autoSave;
        end
        
        function off3_computeReducedSpace(this)
            % Offline phase 3: Generate state space reduction
            time = tic;
            if ~isempty(this.SpaceReducer)
                fprintf('Computing reduced space...\n');
                [V, W] = this.SpaceReducer.generateReducedSpace(this);
                this.Data.V = data.FileMatrix(V,'Dir',this.Data.DataDirectory);
                if ~isempty(W)
                    this.Data.W = data.FileMatrix(W,'Dir',this.Data.DataDirectory);
                else 
                    this.Data.W = this.Data.V;
                end
            else
                this.Data.V = [];
                this.Data.W = [];
            end
            this.OfflinePhaseTimes(3) = toc(time);
            this.autoSave;
        end
        
        function off4_genApproximationTrainData(this)
            % Generates the training data `x_i` for the
            % `\hat{f}`-approximation and precomputes `f(x_i)` values.
            %
            % @todo
            % - include/check MultiArgumentEvaluations-possibility in parallel execution code
            % - Deactivated due to immense overhead and matlab crashes. investigate further
            t = tic;
            if ~isempty(this.Approx)
                this.Data.ApproxTrainData = data.ApproxTrainData.computeFrom(this, ...
                    this.System.f, this.Approx.TrainDataSelector, this.ComputeParallel);
            end
            this.OfflinePhaseTimes(4) = toc(t);
            this.autoSave;
        end
        
        function off5_computeApproximation(this)
            % Offline phase 5: Core function approximation
            %
            % @todo proper subset selection algorithm/method
            time = tic;
            
            if ~isempty(this.Approx)
                if isempty(this.Data.ApproxTrainData.xi)
                    error('No approximation training data available. Called off4_genApproximationTrainData?');
                end
                
                if ~isempty(this.preApproximationTrainingCallback)
                    this.preApproximationTrainingCallback();
                end
                
                if KerMor.App.Verbose > 0
                    fprintf('Starting approximation process for %d dimensions ...\n',size(this.Data.ApproxTrainData.fxi,1));
                end
                
                warning off MATLAB:nearlySingularMatrix;
                %% Approximate core function (is parallelizable for its own)
                % Compile necessary data
                this.Approx.approximateSystemFunction(this);
                warning on MATLAB:nearlySingularMatrix;
                
                if ~isempty(this.postApproximationTrainingCallback)
                    this.postApproximationTrainingCallback();
                end
            end
            
            this.OfflinePhaseTimes(5) = toc(time);
            this.autoSave;
        end
        
        function off6_prepareErrorEstimator(this)
            % Prepares offline data for a possible error estimator.
            t = tic;
            e = this.ErrorEstimator;
            if ~isempty(e)
                if KerMor.App.Verbose > 0
                    fprintf('Starting error estimator offline computations...\n');
                end
                e.offlineComputations(this);
            end
            this.OfflinePhaseTimes(6) = toc(t);
            this.autoSave;
        end
        
        function offlineGenerations(this)
            % Performs all large offline computations for model reduction.
            %
            % This method is mainly used to compile all large data for the
            % reduction process. Its separation from the buildReducedModel
            % method is only for separation reasons.
            % It calls all of the offX_ - methods in their order.
            %
            % See also: buildReducedModel
            this.OfflinePhaseTimes = [];
            
            this.off1_createParamSamples;
            this.off2_genTrainingData;
            this.off3_computeReducedSpace;
            this.off4_genApproximationTrainData;
            this.off5_computeApproximation;
            this.off6_prepareErrorEstimator;
        end
        
        function [reduced,time] = buildReducedModel(this, target_dim)
            % Builds a reduced model from a full model.
            %
            % This is a convenience method, that only calls the
            % models.ReducedModel constructor passing this model instance
            % and any specified target dimension. Also returns the total
            % time that was needed to create the reduced model.
            %
            % Before calling this method ensure that offlineGenerations was
            % called at least once to provide the model's necessary data
            % for the reduction process.
            %
            % Parameters:
            % target_dim: The target dimension `d` of the reduced model.
            % Uses the first `d` columns of the projection matrices `V`
            % (and `W` if set, respectively) to create a subspace-projected
            % reduced model.
            %
            % Return values:
            % reduced: The reduced model created from this full model.
            % @type models.ReducedModel
            % time: The time needed to build the reduced model. @type
            % double
            %
            % See also: offlineGenerations
            tic;
            if nargin < 2
                target_dim = size(this.Data.V,2);
            end
            reduced = models.ReducedModel(this, target_dim);
            time = toc;
        end
        
        function [t, x, time, cache] = computeTrajectory(this, mu, inputidx)
            % Computes a solution/trajectory for the given mu and inputidx in the SCALED state
            % space.
            %
            % Parameters:
            % mu: The parameter `\mu` for the simulation @type colvec<double>
            % inputidx: The integer index of the input function to use. If
            % more than one inputs are specified this is a necessary
            % argument.
            %
            % Return values:
            % t: The times at which the model was evaluated. Will equal
            % the property Times
            % x: The state variables at the corresponding times t.
            % cached: Flag that indicates that this trajectory is loaded from cache. Possible
            % values are 0: Found in no cache, 1: found in ModelData.TrajectoryData, 2: Found
            % in ModelData.SimCache @type integer
            
            if KerMor.App.UseDPCM
                DPCM.criticalsCheck(this);
            end
            if ~isempty(this.System.f) && (isempty(this.System.f.xDim) || isempty(this.System.f.fDim))
                error('Nonlinarity properties fDim and xDim not set. See dscomponents.ACoreFun');
            end
            
            % Try local model data first
            cache = 1;
            [x, time] = this.Data.TrajectoryData.getTrajectory(mu, inputidx);
            if this.EnableTrajectoryCaching && isempty(x)
                cache = 2;
                [x, time] = this.Data.SimCache.getTrajectory([this.T; this.dt; mu], inputidx);
            end
            if ~isempty(x)
                t = this.scaledTimes;
                if size(x,2) ~= length(t)
                    error(['Inconsistent trajectory data! Size of (mu,inputidx)-matching trajectory differs from the model''s Times property.\n'...
                           'Did you change model.dt or model.T and leave model.Data filled with old trajectories?']);
                end
            else
                st = tic;
                [t, x] = computeTrajectory@models.BaseModel(this, mu, inputidx);
                time = toc(st);
                if this.EnableTrajectoryCaching
                    this.Data.SimCache.addTrajectory(x, [this.T; this.dt; mu], inputidx, time);
                end
                cache = 0;
            end
        end
        
        function file = save(this)
            % Saves this instance inside the data.ModelData.DataDirectory folder using the
            % model's SaveTag if set or "model.mat".
            %
            % @new{0,7,dw,2013-05-28} Added this method
            if ~isempty(this.SaveTag)
                filename = this.SaveTag;
            else
                warning('KerMor:NoSaveTag','No SaveTag set. Using default "model.mat"');
                filename = 'model.mat';
            end
            file = fullfile(this.Data.DataDirectory,filename);
            m = this; %#ok (renaming 
            save(file,'m');
        end
    end
    
    methods(Access=private)
        function intTimeChanging(this, src, ~)
            this.ftcold.(src.Name) = this.(src.Name);
        end
        
        function intTimeChanged(this, src, ~)
            if this.ftcold.(src.Name) ~= this.(src.Name)
                md = this.Data.TrajectoryData;
                if md.getNumTrajectories > 0
                    %md.clearTrajectories;
                    %fprintf(2,'Deleted all old trajectories for %s=%f\n',src.Name,this.ftcold.(src.Name));
                end
                %this.Data.SimCache.clearTrajectories;
            end
        end
        
        function nameChanged(this, ~, ~)
            if isempty(this.SaveTag)
                this.SaveTag = regexprep(...
                    strrep(strtrim(this.Name),' ','_'),'[^\d\w~-]','');
            end
        end
        
        function autoSave(this)
            % Only save if in global offline phase
            if this.AutoSave
                save(fullfile(this.Data.DataDirectory,'model_autosave.mat'),'this');
            end
        end
    end
        
    %% Getter & Setter
    methods
        function set.Sampler(this, value)
            this.checkType(value, 'sampling.BaseSampler');
            this.Sampler = value;
        end

        function set.Data(this, value)
            this.checkType(value, 'data.ModelData');
            if ~isequal(this.Data, value)
                this.Data = [];
                this.Data = value;
            end
        end

        function set.SpaceReducer(this, value)
            this.checkType(value, 'spacereduction.BaseSpaceReducer');
            this.SpaceReducer = value;
        end

        function set.Approx(this, value)
            this.checkType(value, 'approx.BaseApprox');
            this.Approx = value;
        end

        function set.preApproximationTrainingCallback(this, value)
            if ~isempty(value) && ~isa(value, 'function_handle')
                error('Value must be a function handle taking the current model instance');
            end
            this.preApproximationTrainingCallback = value;
        end

        function set.postApproximationTrainingCallback(this, value)
            if ~isempty(value) && ~isa(value, 'function_handle')
                error('Value must be a function handle taking the current model instance');
            end
            this.postApproximationTrainingCallback = value;
        end

        function set.TrainingInputs(this, value)
            if ~isempty(value)
                if any(value < 1)
                    error('Value may only contain valid indices for the Inputs cell array.');
                elseif any(value > this.System.InputCount) || any(value < 1)
                    error('Invalid indices for Inputs.');
                elseif isempty(this.System.B)
                    error('You must set the system''s property B when using training inputs.');
                end
            end
            this.fTrainingInputs = value;
        end
        
        function set.EnableTrajectoryCaching(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('EnableTrajectoryCaching must be a flag (scalar logical)');
            end
            this.EnableTrajectoryCaching = value;
        end
        
        function set.AutoSave(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('AutoSave must be a flag (scalar logical)');
            end
            this.AutoSave = value;
        end

        function ti = get.TrainingInputs(this)
            ti = this.fTrainingInputs;
            %             if isempty(ti) && ~isempty(this.System)
            %                 ti = 1:this.System.InputCount;
            %             end
        end

        function c = get.TrainingInputCount(this)
            c = length(this.TrainingInputs);
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this, sobj)
            if ~isa(this,'models.BaseFullModel') && nargin < 2
                error('The model class has changed but the loadobj method does not implement backwards-compatible loading behaviour.\nPlease implement the loadobj-method in your subclass and pass the loaded object struct as second argument.');
            end
            if nargin == 2
                this.Data = sobj.Data;
                this.SpaceReducer = sobj.SpaceReducer;
                this.Approx = sobj.Approx;
                this.JKerMorExport = sobj.JKerMorExport;
                this.Sampler = sobj.Sampler;
                this.EnableTrajectoryCaching = sobj.EnableTrajectoryCaching;
                this.preApproximationTrainingCallback = sobj.preApproximationTrainingCallback;
                this.postApproximationTrainingCallback = sobj.postApproximationTrainingCallback;
                this.ErrorEstimator = sobj.ErrorEstimator;
                this.ComputeTrajectoryFxiData = sobj.ComputeTrajectoryFxiData;
                this.SaveTag = sobj.SaveTag;
                this.fTrainingInputs = sobj.fTrainingInputs;
                this.OfflinePhaseTimes = sobj.OfflinePhaseTimes;
                this.ftcold = sobj.ftcold;
                if isfield(sobj,'AutoSave')
                    this.AutoSave = logical(sobj.AutoSave);
                end
            end
            this = loadobj@models.BaseModel(this);
            this.addlistener('T','PreSet', @this.intTimeChanging);
            this.addlistener('dt','PreSet', @this.intTimeChanging);
            this.addlistener('T','PostSet', @this.intTimeChanged);
            this.addlistener('dt','PostSet', @this.intTimeChanged);
            this.addlistener('Name','PostSet', @this.nameChanged);
        end
    end
           
    methods(Static)
        
        function test_BaseModels
            m = models.BaseFullModel;
            af = dscomponents.AffLinCoreFun;
            af.addMatrix('1', rand(4,4));
            af.addMatrix('sum(mu)*5', rand(4,4)*3);
            
            m.System.f = af;
            m.System.x0 = dscomponents.ConstInitialValue(sin(1:4)');
            
            m.offlineGenerations;
            red = m.buildReducedModel;
            % Test simulations
            m.simulate();
            red.simulate();
            
            % dont forget to free resources! (static handles seem to
            % persist over several method calls)
            clear af m red;
        end
        
        function res = test_BareModel
            m = models.BaseFullModel;
            m.SpaceReducer = [];
            m.Approx = [];
            m.Sampler = [];
            A = rand(2,2);
            m.System = models.BaseDynSystem(m);
            m.System.f = dscomponents.PointerCoreFun(@(x,t,mu)A*x,2);
            m.System.x0 = dscomponents.ConstInitialValue(sin(1:2)');
            m.simulate();
            m.offlineGenerations;
            res = false;
            try 
                m.buildReducedModel;
            catch ME%#ok
                res = true;
            end
            
            % dont forget to free resources! (static handles seem to
            % persist over several method calls)
            clear af m;
        end
        
        function test_LinearModel
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Algorithm.ExpConfig.ParamConfig = [];
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin,ts.testdim);
            s.x0 = ts.x0;
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
        
        function test_LinearModelNoProj
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Algorithm.ExpConfig.ParamConfig = [];
            m.SpaceReducer = [];
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin,ts.testdim);
            s.x0 = ts.x0;
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
        
        function test_LinearModelNoApprox
            ts = testing.testsettings;
            m = ts.m;
            m.Approx = [];
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin,ts.testdim);
            s.x0 = ts.x0;
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
        
        function test_LinearModelParams
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Expansion.ParamKernel = ts.ParamKernel;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin_p,ts.testdim);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            m.simulate(m.getRandomParam);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.getRandomParam);
            clear s m r;
        end
        
        function test_LinearModelInputs
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Algorithm.ExpConfig.ParamConfig = [];
            s = m.System;
            s.B = dscomponents.PointerInputConv(ts.B);
            s.f = dscomponents.PointerCoreFun(ts.flin,ts.testdim);
            s.x0 = ts.x0;
            s.Inputs = ts.Inputs;
            m.simulate([],1);
            m.TrainingInputs = 1:ts.testinputs;
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate([],1);
            clear s m r;
        end
        
        function test_LinearModelParamsInput
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin_p, ts.testdim);
            s.B = dscomponents.PointerInputConv(ts.B_p);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            s.Inputs = ts.Inputs;
            m.simulate(m.getRandomParam, 1);
            m.TrainingInputs = 1:ts.testinputs;
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.getRandomParam, 1);
            clear s m r;
        end
        
        function test_NonlinearModel
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.fnlin,ts.testdim);
            s.x0 = ts.x0;
            m.Approx.Algorithm.ExpConfig.ParamConfig = [];
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
        
        function test_NonlinearModelParams
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.fnlin_p,ts.testdim);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            m.simulate(m.getRandomParam);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.getRandomParam);
            clear s m r;
        end
        
        function test_NonlinearModelInputs
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Algorithm.ExpConfig.ParamConfig = [];
            s = m.System;
            s.B = dscomponents.PointerInputConv(ts.B);
            s.f = dscomponents.PointerCoreFun(ts.fnlin,ts.testdim);
            s.x0 = ts.x0;
            s.Inputs = ts.Inputs;
            m.simulate([],1);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate([],1);
            clear s m r;
        end
        
        function test_NonlinearModelParamsInput
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Expansion.ParamKernel = ts.ParamKernel;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.fnlin_p,ts.testdim);
            s.B = dscomponents.PointerInputConv(ts.B_p);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            s.Inputs = ts.Inputs;
            m.simulate(m.getRandomParam,1);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.getRandomParam,1);
            clear s m r;
        end
        
        function test_TimeDependentOutput
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.Algorithm.ExpConfig.ParamConfig = [];
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin,ts.testdim);
            s.x0 = ts.x0;
            s.C = dscomponents.PointerOutputConv(@(t,mu)t,true);
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
        
        function res = test_FileSystemTidyness
            % Plain small model
            m = models.circ.RCLadder(10);
            m.offlineGenerations;
            dir = m.Data.DataDirectory;
            
            % Plain small model with FileData
            m = models.circ.RCLadder(10);
            
            % Overwrite test: old m instance isnt needed anywhere else
            res = exist(dir,'file') == 0;
            
            % Test with FileTrajectoryData, FileTrajectoryFxiData and nonempty SimCache
            m.Data.useFileTrajectoryData;
            m.ComputeTrajectoryFxiData = true;
            m.TrainingInputs = [2 3];
            m.offlineGenerations;
            m.simulate([],4);
            
            % "Tweak" atd FileMatrices to contain more than one block
            atd = m.Data.ApproxTrainData;
            atd.xi = atd.xi.copyWithNewBlockSize(0.1);
            atd.fxi = atd.fxi.copyWithNewBlockSize(0.1);
            
            dir = m.Data.DataDirectory;
            clear atd;
            clear m;
            res = res && exist(dir,'file') == 0;
        end
    end
end
