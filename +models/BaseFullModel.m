classdef BaseFullModel < models.BaseModel & IParallelizable
% The base class for any KerMor detailed model
%
% Implementers of custom models are to inherit from this base class
% in order for it to work with KerMor.
% For custom models, the properties of this class (combined with
% those from BaseModel) can be set to influence the model behaviour
% and reduction methods.
% For the implementation of custom dynamical systems, refer to
% BaseFirstOrderSystem.
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
% See also: models BaseModel BaseFirstOrderSystem
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        
        % Flag to determine if the data collected by the
        % approx.Approx.TrainDataSelector should be projected into the
        % previously computed projection space V (see offline phase 3).
        % If so, the `f(x_i)` values are also re-computed for the updated
        % `x_i` samples.
        %
        % This feature might be of importance when kernel methods are used
        % for nonlinearity approximation, and the projection requires the
        % centers to reside in the projection space V.
        %
        % As the more common choice in KerMor so far is DEIM, the default
        % is false as `x_i` data is not needed for DEIM.
        %
        % See also: data.ApproxTrainData approx.DEIM
        %
        % @type logical @default false
        ProjectApproxTrainingData = false;
        
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
        % The indices of inputs to use for training data generation. Uses
        % the DefaultInput if not set (and DefaultInput is not empty)
        %
        % @propclass{optional} Set subindices to restrict offline generation phase to a subset of
        % inputs.
        %
        % @default DefaultInput @type integer
        TrainingInputs;
        
        % The indices of the model parameters to use for training data
        % generation.
        %
        % The model's DefaultMu will be used as "base" parameter, within
        % which the TrainingParams will be changed accordingly.
        %
        % @propclass{optional} Set subindices to restrict offline generation phase to a subset of
        % parameters.
        %
        % @default 1:ParamCount @type rowvec<integer>
        TrainingParams;
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
        
        % Contains a flag for each input/parameter combination which has
        % been successfully computed completely during the offline phase.
        %
        % This output will depend on the UniformTrajectories setting of the
        % ModelData.TrajectoryData instances; if set to true, offline phase
        % 2 will not run through anyways.
        TrajectoriesCompleted = [];
    end
    
    properties(Access=private)
        fTrainingInputs = [];
        fTrainingParams = [];
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
            this.System = models.BaseFirstOrderSystem(this);
            
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
            time = tic;
            
            num_s = max(1,this.Data.SampleCount);
            num_in = max(1,this.TrainingInputCount);
            
            completed = zeros(3,num_s*num_in);
            tlen = length(this.Times);
            
            % Clear old trajectory data.
            this.Data.TrajectoryData.clearTrajectories;
            this.RealTimePlotting = false;
            
            oldv = this.EnableTrajectoryCaching;
            this.EnableTrajectoryCaching = false;
            
            % Automatically use the default input if exists and no others
            % are set.
            if isempty(this.TrainingInputs) && ~isempty(this.DefaultInput)
                this.TrainingInputs = this.DefaultInput;
            end
            
            %% Parallel - computation
            if this.ComputeParallel
                told = [];
                if ~isa(this.Data.TrajectoryData,'data.FileTrajectoryData')
                    told = this.Data.TrajectoryData;
                    this.Data.useFileTrajectoryData;
                end
                idxmat = Utils.createCombinations(1:num_s,this.TrainingInputs);
                
                fprintf('Starting parallel projection training data computation of %d trajectories on %d workers...\n',...
                    size(idxmat,2),PCPool.size);
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
                    [t, x, ctime] = remote.computeTrajectory(mu, inputidx);
                    
                    % Assign snapshot values
                    remote.Data.TrajectoryData.addTrajectory(x, mu, inputidx, ctime);

                    % Evaluate fxi on current values if required
                    if this.ComputeTrajectoryFxiData
                        hlp = tic;
                        fx = this.System.f.evaluateMulti(x, t, mu);
                        ctime = toc(hlp);
                        this.Data.TrajectoryFxiData.addTrajectory(fx, mu, inputidx, ctime);
                    end
                    completed(:,idx) = [inidx(idx); paridx(idx); size(x,2) == tlen];
                end
                % Build the hash map for the local FileTrajectoryData (parfor
                % add them to remotely instantiated FileTrajectoryDatas and does
                % not sync them (not generically possible)
                
                td = this.Data.TrajectoryData;
                td.consolidate(this);
                if ~isempty(told)
                    told.transferFrom(td);
                    td.clearTrajectories;
                    this.Data.TrajectoryData = [];
                    this.Data.TrajectoryData = told;
                end
                
                if ~isempty(this.Data.TrajectoryFxiData) && isa(this.Data.TrajectoryFxiData,'data.FileTrajectoryData')
                    this.Data.TrajectoryFxiData.consolidate(this);
                end
                
                %% Non-parallel computation
            else
                % Assume no parameters or inputs
                mu = [];
                inputidx = [];
                
                if KerMor.App.Verbose > 0
                    pi = ProcessIndicator(sprintf('Generating projection training data (%d trajectories)...',num_in*num_s),num_in*num_s);
                end
                
                completed = double.empty(3,0);
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
                        
                        % Get trajectory
                        try
                            [t, x, ctime, cache] = this.computeTrajectory(mu, inputidx);
                            % Store snapshot values, only if not already from trajectory data
                            if cache ~= 1
                                if size(x,2) == tlen
                                    this.Data.TrajectoryData.addTrajectory(x, mu, inputidx, ctime);
                                elseif size(x,2) ~= tlen
                                    this.Data.TrajectoryIncompleteData.addTrajectory(x, mu, inputidx, ctime);
                                end
                            end

                            % Evaluate fxi on current values if required
                            if this.ComputeTrajectoryFxiData
                                hlp = tic;
                                fx = this.System.f.evaluateMulti(x, t, mu);
                                ctime = toc(hlp);
                                this.Data.TrajectoryFxiData.addTrajectory(fx, mu, inputidx, ctime);
                            end
                            ok = size(x,2) == tlen;    
                        catch ME
                            ok = 0;
                            warning(ME.identifier,'Failed to compute trajectory:\n%s',ME.message);
                        end
                        completed(:,end+1) = [inidx; pidx; ok]; %#ok
                        
                        if KerMor.App.Verbose > 0
                            pi.step;
                        end
                    end
                end
                if KerMor.App.Verbose > 0
                    pi.stop;
                end
            end
            
            this.EnableTrajectoryCaching = oldv;
            this.TrajectoriesCompleted = completed;
            
            % delete parameters of incomplete trajectories
            incomplete = find(completed(3,:)==0);
            this.Data.ParamSamples(:,incomplete) = [];
            
            % Input space span computation (if spacereduction is used)
            if ~isempty(this.SpaceReducer) && any([this.SpaceReducer(:).IncludeBSpan])
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
            this.Data.ProjectionSpaces = [];
            if ~isempty(this.SpaceReducer)
                ns = length(this.SpaceReducer);
                %% Sanity checks
                dim = this.System.NumTotalDofs;
                numalgdofs = this.System.NumAlgebraicDofs;
                for k = 1:ns
                    td = this.SpaceReducer(k).TargetDimensions;
                    if strcmp(td,':')
                        if numalgdofs > 0
                            this.SpaceReducer(k).TargetDimensions = 1:dim-numalgdofs;
                            fprintf('Adjusting TargetDimensions of space reducer#%d: Selection '':'' invalid when having algebraic constraints.\n',k);
                        end
                    elseif any(td > dim-numalgdofs)
                        disp(td(td > dim-numalgdofs));
                        error('Cannot select algebraic constraints for reduction (Check System.AlgebraicConditionDoF and SpaceReducer.TargetDimensions) ! Conflicting DoF indices above.');
                    end
                end
                %% Actual computation
                fprintf('Computing %d projection space(s)...\n',ns);
                for k = 1:ns
                    s = this.SpaceReducer(k);
                    [V, W] = s.generateReducedSpace(this);
                    this.Data.addProjectionSpace(V,W,s.TargetDimensions);
                end
            end
            this.OfflinePhaseTimes(3) = toc(time);
            this.autoSave;
        end
        
        function off4_genApproximationTrainData(this)
            % Generates the training data `x_i` for the
            % `\hat{f}`-approximation and precomputes `f(x_i)` values.
            %
            % The main functionality of this phase is located at the
            % ApproxTrainData class as it is used more than once and not
            % exclusively at this location
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
                
                %warning off MATLAB:nearlySingularMatrix;
                %% Approximate core function (is parallelizable for its own)
                % Compile necessary data
                this.Approx.approximateSystemFunction(this);
                %warning on MATLAB:nearlySingularMatrix;
                
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
        
        function [reduced,time] = buildReducedModel(this, varargin)
            % Builds a reduced model from a full model.
            %
            % Before calling this method ensure that offlineGenerations was
            % called at least once to provide the model's necessary data
            % for the reduction process.
            %
            % Parameters:
            % target_dim: The target dimensions `d` of the reduced model.
            % If more than one subspace is set to be computed, target_dim
            % may be a vector determining the desired size of each
            % subspace. The overall dimension will then be the sum over all
            % single ones.
            % If algebraic DoFs are present, the reduced space will always
            % have all of those added to it's dimensions
            %
            % Return values:
            % reduced: The reduced model created from this full model.
            % @type models.ReducedModel
            % time: The time needed to build the reduced model. @type
            % double
            %
            % See also: offlineGenerations
            tic;
            reduced = models.ReducedModel(this);
            reduced.build(varargin{:});
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
            if isempty(x) && this.EnableTrajectoryCaching
                cache = 2;
                [x, time] = this.Data.SimCache.getTrajectory([this.T; this.dt; mu], inputidx);
            end
            if ~isempty(x)
                t = this.scaledTimes;
                if size(x,2) < length(t)
                    warning('Inconsistent trajectory data: Cached trajectory shorter than current time-steps. Using matching subarray of model.scaledTimes');
                   t = t(1:size(x,2));   
                elseif size(x,2) > length(t)
                    error('Inconsistent trajectory data: More trajectory samples than time-steps in current model. Please fix.');
                end
                % Also set the current mu config for output mapping &
                % plotting to work correctly.
                this.System.setConfig(mu, inputidx);
            else
                st = tic;
                
                % Preparation call to approximation, if present
                if ~isempty(this.Approx)
                    this.Approx.prepareSimulation(mu);
                end
                
                [t, x] = computeTrajectory@models.BaseModel(this, mu, inputidx);
                
                time = toc(st);
                if this.EnableTrajectoryCaching
                    if size(x,2) == length(this.Times)
                        this.Data.SimCache.addTrajectory(x, [this.T; this.dt; mu], inputidx, time);
                    else
                        warning('KerMor:InconsistentTrajectoryData',...
                            'Not caching trajectory: Returned result size differs from simulation times.');
                    end
                end
                cache = 0;
            end
        end
        
        function [t, y, sec, lastx] = simulateIntervals(this, mu, inputidx, chunksize, callback)
            % Simulate the current model in intervals and return only a
            % selected output.
            %
            % This convenience method allows to perform simulations where
            % the full state space vector would not fit into memory.
            %
            % Parameters:
            % mu: The desired model parameter @type colvec<double>
            % inputidx: The input to use. Set to [] for no input. @type
            % handle
            % chunksize: The maximum size of a full single trajectory in
            % GB @type double @default Quarter of the system's memory
            % callback: A function handle to use for extracting the
            % quantities of interest from each interval solution.
            % Caution: There are no checks as to whether the callback tries
            % to create matrices that are too large for the current system.
            % Usually, the callback extracts the system output or a
            % small-dimensional subset @type function_handle
            % @default this.System.computeOutput
            %
            % Return values:
            % t: The computed total times
            % y: The total output, where each interval's result has been
            % processed by the callback to obtain y.
            % sec: The total computation time in seconds @type double
            % lastx: The state of the full system at the final time T. Can
            % be used as e.g. initial value for further simulations.
            % @type colvec<double>
            
            if nargin < 5
                callback = @this.System.computeOutput;
                if nargin < 4
                    mem = memory;
                    % Default to a quarter of the available RAM on this machine
                    chunksize = mem.TotalSystemMemory/1024^3/4;
                end
            end
            t = this.Times;
            total_len = length(t);
            nchunks = ceil(this.FullStateTrajectorySize/chunksize);
            if nchunks == 1
                [~,~,sec,x] = this.simulate(mu,inputidx);
                hlp = tic;
                y = callback(x);
                sec = sec+toc(hlp);
            else
                interval_idx = floor(linspace(1,total_len,nchunks+1));
                intervals = this.Times(interval_idx);
                numint = length(intervals)-1;
                % Store original values
                oldT = this.T;
                oldt0 = this.t0;
                oldx0 = this.System.x0;
                tc = this.EnableTrajectoryCaching;
                this.EnableTrajectoryCaching = false;
                try
                    sec = 0;
                    for iidx = 1:numint
                        this.t0 = intervals(iidx);
                        this.T = intervals(iidx+1);
                        fprintf('Starting interval %d of %d (size %gMB) ... ',...
                            iidx,numint,this.FullStateTrajectorySize*1024);
                        [~,~,ct,x] = this.simulate(mu,inputidx);
                        hlp = tic;
                        y(:,interval_idx(iidx):interval_idx(iidx+1)) = callback(x);
                        sec = sec+ct+toc(hlp);
                        this.System.x0 = dscomponents.ConstInitialValue(x(:,end));
                    end
                    lastx = x(:,end);
                catch ME
                    this.T = oldT;
                    this.t0 = oldt0;
                    this.System.x0 = oldx0;
                    this.EnableTrajectoryCaching = tc;
                    rethrow(ME);
                end
                this.T = oldT;
                this.t0 = oldt0;
                this.System.x0 = oldx0;
                this.EnableTrajectoryCaching = tc;
            end
        end
        
        function [t,y,mu,in,ct] = getSimCacheTrajectory(this, nr)
            [y,mu,in,ct] = this.Data.SimCache.getTrajectoryNr(nr);
            this.T = mu(1);
            this.dt = mu(2);
            mu = mu(3:end);
            t = this.Times;
            this.System.setConfig(mu,in);
        end
        
        function plotTrajectoryNr(this, nr)
            [y,mu,in] = this.Data.TrajectoryData.getTrajectoryNr(nr);
            this.System.setConfig(mu,in);
            this.plot(this.Times,y);
        end
        
        function [V, W] = assembleProjectionMatrices(this, target_dim)
            fd = this.Data;
            if isempty(fd.ProjectionSpaces)
                V = [];
            else
                if nargin < 2
                    % Default: create a reduced model 10th the size of the full
                    % one
                    target_dim = max(1,floor(this.System.NumTotalDofs/10));
                end
                if target_dim > size(fd.ProjectionSpaces(1).V,2)
                    warning('Target dimension %d too large, using max value of %d.',...
                        target_dim,size(fd.ProjectionSpaces(1).V,2));
                    target_dim = size(fd.ProjectionSpaces(1).V,2);
                end
                V = fd.ProjectionSpaces(1).V(:,1:target_dim);
            end
            % We do 
            W = V;
            %% Old code for multiple projection spaces (not currently used anywhere)
%             md = this.Data;
%             ps = md.ProjectionSpaces;
%             ns = length(ps);
%             
%             % If no target dimension(s) are given, assume full size (for
%             % each subspace)
%             if nargin < 2 || isempty(target_dim)
%                 for k = 1:ns
%                     target_dim(k) = ps(k).Size;
%                 end
%                 % If only a scalar is given, assume a total desired size of
%                 % target_dim and equally split sizes for each subspace
%             elseif isscalar(target_dim)
%                 target_dim = ones(ns,1)*floor(target_dim/ns);
%             end
%             for k = 1:ns
%                 if target_dim(k) > ps(k).Size
%                     warning('Target size %d for subspace %d too large! Using max value of %d',...
%                         target_dim(k),k,ps(k).Size);
%                     target_dim(k) = ps(k).Size;
%                 end
%                 % Store the currently used effective size of that subspace
%                 % for use by other components upon projection
%                 ps(k).LastEffectiveSize = target_dim(k);
%             end
%             
%             %% Compute the effective projection matrix
%             % As we might have multiple projection subspaces, we need to
%             % assemble a global projection matrix (block-like if separate
%             % subspace dimensions are "sorted" in full model).
%             % Any not-projected DoFs are forwarded as-is by identity.
%             % DAE constraints are treated outside of this method - this is
%             % purely dof-related.
%             V = []; W = []; %#ok<*PROP>
%             if ns > 0
%                 dim = this.System.NumTotalDofs;
%                 V = zeros(dim,sum(target_dim));
%                 W = V;
%                 offset = 0;
%                 done = false(dim,1);
%                 for k=1:ns
%                     s = ps(k);
%                     sel = 1:target_dim(k);
%                     V(s.Dimensions,offset + sel) = s.V(:,sel);
%                     if ~isempty(s.W)
%                         W(s.Dimensions,offset + sel) = s.W(:,sel);
%                     end
%                     offset = offset + target_dim(k);
%                     done(s.Dimensions) = true;
%                 end
% %                 % Insert identity for any remaining dimensions (corresponds
% %                 % to algebraic constraint DoFs
% %                 unreduced = find(~done);
% %                 nunred = length(unreduced);
% %                 V(unreduced,end+(1:nunred)) = eye(nunred);
% %                 W(unreduced,end+(1:nunred)) = eye(nunred);
%             end
        end
        
        function file = save(this, directory, filename)
            % Saves this instance inside the data.ModelData.DataDirectory folder using the
            % model's SaveTag if set or "model.mat".
            %
            % Parameters:
            % directory: The directory to save to. @type char @default
            % data.ModelData.DataDirectory
            % filename: The filename to use. Will derive from
            % models.BaseFullModel.SaveTag if not specified, and if the
            % SaveTag is not set, 'model.mat' will be used. @type char
            % @default models.BaseFullModel.SaveTag
            %
            % @new{0,7,dw,2013-05-28} Added this method
            if nargin < 3
                if ~isempty(this.SaveTag)
                    filename = this.SaveTag;
                else
                    warning('KerMor:NoSaveTag','No SaveTag set. Using default "model.mat"');
                    filename = 'model.mat';
                end
            end
            if nargin < 2
                directory = this.Data.DataDirectory;
            end
            file = fullfile(directory,filename);
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
        
        function set.TrainingParams(this, value)
            if ~isempty(value)
                if any(value < 1)
                    error('Value may only contain valid indices for the parameter array.');
                elseif any(value > this.System.ParamCount) || any(value < 1)
                    error('Invalid indices for Inputs.');
                end
            end
            this.fTrainingParams = value;
        end
        
        function set.ProjectApproxTrainingData(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('ProjectApproxTrainingData must be a flag (scalar logical)');
            end
            this.ProjectApproxTrainingData = value;
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
        end
        
        function ti = get.TrainingParams(this)
            ti = this.fTrainingParams;
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
                this = loadobj@models.BaseModel(this, sobj);
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
                if isfield(sobj,'TrajectoriesCompleted')
                    this.TrajectoriesCompleted = sobj.TrajectoriesCompleted;
                end
                % Load data class at last; any other components that might
                % make use of stuff inside the model folder should be able
                % to clean up
                this.Data = sobj.Data;
            else
                this = loadobj@models.BaseModel(this);
            end
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
            af = dscomponents.AffLinCoreFun(m.System);
            af.addMatrix('1', rand(4,4));
            af.addMatrix('sum(mu)*5', rand(4,4)*3);
            
            m.System.f = af;
            m.System.x0 = dscomponents.ConstInitialValue(sin(1:4)');
            m.System.addParam('test',.5,'Range',[0 1],'Desired',10);
            m.System.MaxTimestep = m.dt/2;
            
            m.offlineGenerations;
            red = m.buildReducedModel;
            % Test simulations
            m.simulate(m.getRandomParam);
            red.simulate(m.getRandomParam);
            
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
            m.System = models.BaseFirstOrderSystem(m);
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
                s.addParam(p.Name, p.Default, 'Range', [p.MinVal p.MaxVal], 'Desired', p.Desired);
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
                s.addParam(p.Name, p.Default, 'Range', [p.MinVal p.MaxVal], 'Desired', p.Desired);
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
                s.addParam(p.Name, p.Default, 'Range', [p.MinVal p.MaxVal], 'Desired', p.Desired);
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
                s.addParam(p.Name, p.Default, 'Range', [p.MinVal p.MaxVal], 'Desired', p.Desired);
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
