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
% See also: models BaseModel BaseDynSystem
%
% @author Daniel Wirtz @date 16.03.2010
%
% @change{0,6,dw,2012-04-26} Using tools.ProcessIndicator for training data
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
% data.FileModelData (parallel execution did not sync the hashmaps, now
% running data.FileModelData.consolidate fixes this)
%
% @change{0,5,dw,2011-08-04} Adopted the off2_genTrainingData method to the new data.AModelData
% structure. Now all trajectories are stored either in memory or disk, and the data.AModelData
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
    
    properties
        % The full model's data container.
        % Defaults to an empty container.
        %
        % @propclass{data} The model's data container
        %
        % @default data.MemoryModelData @type data.AModelData
        %
        % See also: data.AModelData
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
        % @default approx.KernelApprox
        % @type approx.BaseApprox
        %
        % See also: approx approx.algorithms approx.selection
        Approx;
        
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
        % @propclass{optional} Caching can speedup computations but may require extra disk
        % space.
        %
        % @type logical @default true
        EnableTrajectoryCaching = true;
        
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
        % Contains a `1\times 5` row vector with the corresponding times for the 'off1_' to 'off_5'
        % phases.
        %
        % @default []
        % @type rowvec
        OfflinePhaseTimes = [];
    end
    
    properties(Access=private)
        fTrainingInputs = [];
        ftcold = struct;
    end
    
    properties(SetAccess=private, GetAccess=protected, Transient)
        simCache;
    end
           
    methods
        
        function this = BaseFullModel
            % Creates a new instance of a full model.
            
            this = this@models.BaseModel;
            
            % Setup default values for the full model's components
            this.Sampler = sampling.RandomSampler;
            
            this.SpaceReducer = spacereduction.PODGreedy;
            
            this.Approx = approx.KernelApprox;
            
            % ModelData defaults to MemoryData
            this.Data = data.MemoryModelData;
            
            % Set default dynamical system
            this.System = models.BaseDynSystem(this);
            
            % Register default properties
            this.registerProps('Sampler','SpaceReducer','Approx',...
                'preApproximationTrainingCallback','postApproximationTrainingCallback',...
                'Data','TrainingInputs');
            
            % Create a listener for the 
            this.addlistener('T','PreSet', @this.intTimeChanging);
            this.addlistener('dt','PreSet', @this.intTimeChanging);
            this.addlistener('T','PostSet', @this.intTimeChanged);
            this.addlistener('dt','PostSet', @this.intTimeChanged);
            this.simCache = data.FileModelData(this, KerMor.App.TempDirectory);
        end
        
        function delete(this)
            % @todo in AModelData: delete datadir if empty on destruction
            % Clean up the simulation cache
            
            % remove all temp files
            %this.simCache.clearTrajectories;
            % remove folder
            %this.simCache.delete;
        end
        
        function time = off1_createParamSamples(this)
            % Offline phase 1: Sample generation.
            
            % Sampling
            time = tic;
            if ~isempty(this.Sampler)
                this.Data.ParamSamples = this.Sampler.generateSamples(this);
            elseif this.System.ParamCount > 0
                error('No parameter sampler set for a parameterized system. See package sampling and the Sampler property of the BaseFullModel class.');
            end
            time = toc(time);
        end
        
        function time = off2_genTrainingData(this)
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
            %                 % @TODO: switch back if changed!
            %                 this.ComputeParallel = false;
            %             end
            
            % Clear old trajectory data.
            this.Data.clearTrajectories;
            this.RealTimePlotting = false;
            
            %% Parallel - computation
            if this.ComputeParallel
                %error('Parallel computing not yet tested with new ModelData structure / model.Data.addTrajectory might not be thread-safe!');
                
                idxmat = general.Utils.createCombinations(1:num_s,1:num_in);
                
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
                    if remote.Data.SampleCount > 0 %#ok<PFBNS>
                        mu = remote.Data.getParams(paridx(idx));
                    end
                    % Check for inputs
                    inputidx = [];
                    if remote.TrainingInputCount > 0
                        inputidx = inidx(idx);
                    end
                    
                    % Get trajectory
                    [~, x] = remote.computeTrajectory(mu, inputidx);
                    
                    % Assign snapshot values
                    remote.Data.addTrajectory(x, mu, inputidx);
                end
                % Build the hash map for the local FileModelData (parfor
                % add them to remotely instantiated FileModelDatas and does
                % not syn them (not generically possible)
                if isa(this.Data,'data.FileModelData')
                    this.Data.consolidate(this);
                end
                
                %% Non-parallel computation
            else
                % Assume no parameters or inputs
                mu = [];
                inputidx = [];
                
                if KerMor.App.Verbose > 0
                    pi = tools.ProcessIndicator(sprintf('Generating projection training data (%d trajectories)...',num_in*num_s),num_in*num_s);
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
                        
                        % Get trajectory (disable caching as trajectories go to this.Data
                        % anyways)
                        oldv = this.EnableTrajectoryCaching;
                        this.EnableTrajectoryCaching = false;
                        [~, x, ctime] = this.computeTrajectory(mu, inputidx);
                        this.EnableTrajectoryCaching = oldv;
                        
                        % Assign snapshot values
                        this.Data.addTrajectory(x, mu, inputidx, ctime);
                        if KerMor.App.Verbose > 0
                            pi.step;
                        end
                    end
                end
                if KerMor.App.Verbose > 0
                    pi.stop;
                end
            end
            time = toc(time);
        end
        
        function time = off3_computeReducedSpace(this)
            % Offline phase 3: Generate state space reduction
            time = tic;
            V = []; W = [];
            if ~isempty(this.SpaceReducer)
                fprintf('Computing reduced space...\n');
                [V, W] = this.SpaceReducer.generateReducedSpace(this);
            end
            this.Data.V = V;
            this.Data.W = W;
            time = toc(time);
        end
        
        function time = off4_genApproximationTrainData(this)
            % Generates the training data `x_i` for the
            % `\hat{f}`-approximation and precomputes `f(x_i)` values.
            %
            % @todo
            % - include/check MultiArgumentEvaluations-possibility in parallel execution code
            % - Deactivated due to immense overhead and matlab crashes. investigate further
            t = tic;
            if ~isempty(this.Approx)
                [atd, time] = computeTrainingData(this,...
                    this.System.f, this.Approx.TrainDataSelector);
                this.Data.ApproxTrainData = atd;
            end
            time = time + toc(t);
        end
        
        function time = off5_computeApproximation(this)
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
            
            time = toc(time);
        end
        
        function time = off6_prepareErrorEstimator(this)
            % Prepares offline data for a possible error estimator.
            t = tic;
            e = this.ErrorEstimator;
            time = 0;
            if ~isempty(e)
                if KerMor.App.Verbose > 0
                    fprintf('Starting error estimator offline computations...\n');
                end
                % Precompute large offline data needed for the DEIM error
                % estimator
                if isa(e,'error.DEIMEstimator')
                    [jtd, time] = computeTrainingData(this,...
                        general.JacCompEvalWrapper(this.System.f),...
                        e.TrainDataSelector);
                    this.Data.JacobianTrainData = jtd;
                    
                    d = this.System.f.XDim;
                    n = size(jtd.fxi,2);
                    v = zeros(d,n);
                    ln = zeros(1,n);
                    times = ln;
                    pi = tools.ProcessIndicator('Computing Jacobian similarity transform data for %d jacobians',n,false,n);
                    for nr = 1:n
                        J = reshape(jtd.fxi(:,nr),d,d);
                        t = tic;
                        [ln(nr), v(:,nr)] = general.Utils.logNorm(J);
                        times(nr) = toc(t);
                        pi.step;
                    end
                    pi.stop;
                    
                    jstd.VFull = v;
                    jstd.LogNorms = ln;
                    jstd.CompTimes = times;
                    this.Data.JacSimTransData = jstd;
                end
                e.offlineComputations(this);
            end
            time = time + toc(t);
        end
        
        function times = offlineGenerations(this)
            % Performs all large offline computations for model reduction.
            %
            % This method is mainly used to compile all large data for the
            % reduction process. Its separation from the buildReducedModel
            % method is only for separation reasons.
            % It calls all of the offX_ - methods in their order.
            %
            % See also: buildReducedModel
            
            times(1) = this.off1_createParamSamples;
            times(2) = this.off2_genTrainingData;
            times(3) = this.off3_computeReducedSpace;
            times(4) = this.off4_genApproximationTrainData;
            times(5) = this.off5_computeApproximation;
            times(6) = this.off6_prepareErrorEstimator;
            
            % Store the computation times
            this.OfflinePhaseTimes = times;
        end
        
        function [reduced,time] = buildReducedModel(this)
            % Builds a reduced model from a full model.
            %
            % Before calling this method ensure that offlineGenerations was
            % called at least once to provide the model's necessary data
            % for the reduction process.
            %
            % Return values:
            % reduced: The reduced model created from this full model.
            % time: The time needed to build the reduced model.
            %
            % See also: offlineGenerations
            % @docupdate
            
            if isempty(this.Data)
                error('No ModelData class found. Forgot to call offlineGenerations?');
            end
            if isempty(this.Data.getNumTrajectories == 0) && ~(this.Data.SampleCount == 0 && this.TrainingInputCount == 0)
                error('No Snapshot data available. Forgot to call offlineGenerations before?');
            end
            tic;
            reduced = models.ReducedModel(this);
            time = toc;
        end
        
        function [t, x, time] = computeTrajectory(this, mu, inputidx)
            % Computes a solution/trajectory for the given mu and inputidx in the SCALED state
            % space.
            %
            % Parameters:
            % mu: The parameter `\mu` for the simulation
            % inputidx: The integer index of the input function to use. If
            % more than one inputs are specified this is a necessary
            % argument.
            %
            % Return values:
            % t: The times at which the model was evaluated. Will equal
            % the property Times
            % x: The state variables at the corresponding times t.
            %
            % @todo check in model data if trajectory is not already present (if FileModelData
            % is used the trajectories of the offline phase are stored twice in m.Data and
            % m.simCache)
            
            if ~isempty(mu) && size(mu,2) > 1
                if size(mu,1) > 1
                    error('The mu parameter must be a single column vector.');
                else
                    warning('KerMor:BaseDynSystem','Please use column vectors for parameters. Reshaping.');
                    mu = reshape(mu,[],1);
                end
            end
            
            if KerMor.App.UseDPCM
                DPCM.criticalsCheck(this);
            end
            
            % Try local model data first
            [x, time] = this.Data.getTrajectory(mu, inputidx);
            if isempty(x) && this.EnableTrajectoryCaching
                [x, time] = this.simCache.getTrajectory([this.T; this.dt; mu], inputidx);
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
                    this.simCache.addTrajectory(x, [this.T; this.dt; mu], inputidx, time);
                end
            end
        end
    end
    
    methods(Access=private)
        function intTimeChanging(this, src, ~)
            this.ftcold.(src.Name) = this.(src.Name);
        end
        
        function intTimeChanged(this, src, ~)
            if ~isempty(this.Data) && this.ftcold.(src.Name) ~= this.(src.Name)
                md = this.Data;
                if md.getNumTrajectories > 0
                    md.clearTrajectories;
                    fprintf(2,'Deleted all old trajectories for %s=%f\n',src.Name,this.ftcold.(src.Name));
                end
            end
        end
        
        function [atd, time] = computeTrainingData(this, f, selector)
            % Internal method that computes training data using a selector
            % and a function
            time = tic;
            
            % Select subset of projection training data
            atd = selector.selectTrainingData(this);

            % If projection is used, train approximating function in
            % centers projected into the subspace.
            if ~isempty(this.Data.V) && ~isempty(this.Data.W)
                atd.xi = this.Data.V*(this.Data.W'*atd.xi);
            end

            % Compute f-Values at training data
            if isempty(atd.mui)
                atdmui = double.empty(0,length(atd.ti));
            else
                atdmui = atd.mui;
            end

            if this.ComputeParallel
                atdxi = atd.xi;
                atdti = atd.ti;
                fval = zeros(size(atdxi));
                if KerMor.App.Verbose > 0
                    fprintf('Starting parallel f-values computation at %d points on %d workers...\n',size(atd,2),matlabpool('size'));
                end
                parfor sidx=1:size(atdxi,2)
                    fval(:,sidx) = ...
                        f.evaluateCoreFun(atdxi(:,sidx),... % x
                        atdti(sidx),... % t
                        atdmui(:,sidx)); %#ok mu
                end
                atd.fxi = fval;
            else
                if KerMor.App.Verbose > 0
                    fprintf('Serial computation of f-values at %d points ...\n',size(atd.xi,2));
                end
                atd.fxi = f.evaluate(atd.xi, atd.ti, atdmui);
            end
            time = toc(time);
        end
    end
        
    %% Getter & Setter
    methods
        function set.Sampler(this, value)
            this.checkType(value, 'sampling.BaseSampler');
            this.Sampler = value;
        end

        function set.Data(this, value)
            this.checkType(value, 'data.AModelData');
            this.Data = value;
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
        function this = loadobj(this)
            this = loadobj@DPCMObject(this);
            this.addlistener('T','PreSet', @this.intTimeChanging);
            this.addlistener('dt','PreSet', @this.intTimeChanging);
            this.addlistener('T','PostSet', @this.intTimeChanged);
            this.addlistener('dt','PostSet', @this.intTimeChanged);
            % Reassign a simulation cache (transient variable)
            this.simCache = data.FileModelData(this, KerMor.App.TempDirectory);
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
            m.offlineGenerations;
            red = m.buildReducedModel;
            red.simulate();
            
            % dont forget to free resources! (static handles seem to
            % persist over several method calls)
            clear af m red;
            res = true;
        end
        
        function test_LinearModel
            ts = testing.testsettings;
            m = ts.m;
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
            m.Approx.ParamKernel = ts.ParamKernel;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin_p,ts.testdim);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            m.simulate(m.System.getRandomParam);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.System.getRandomParam);
            clear s m r;
        end
        
        function test_LinearModelInputs
            ts = testing.testsettings;
            m = ts.m;
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
            m.Approx.ParamKernel = ts.ParamKernel;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin_p, ts.testdim);
            s.B = dscomponents.PointerInputConv(ts.B_p);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            s.Inputs = ts.Inputs;
            m.simulate(m.System.getRandomParam, 1);
            m.TrainingInputs = 1:ts.testinputs;
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.System.getRandomParam, 1);
            clear s m r;
        end
        
        function test_NonlinearModel
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.fnlin);
            s.x0 = ts.x0;
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
        
        function test_NonlinearModelParams
            ts = testing.testsettings;
            m = ts.m;
            m.Approx.ParamKernel = ts.ParamKernel;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.fnlin_p);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            m.simulate(m.System.getRandomParam);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.System.getRandomParam);
            clear s m r;
        end
        
        function test_NonlinearModelInputs
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.B = dscomponents.PointerInputConv(ts.B);
            s.f = dscomponents.PointerCoreFun(ts.fnlin);
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
            m.Approx.ParamKernel = ts.ParamKernel;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.fnlin_p);
            s.B = dscomponents.PointerInputConv(ts.B_p);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            s.Inputs = ts.Inputs;
            m.simulate(m.System.getRandomParam,1);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate(r.System.getRandomParam,1);
            clear s m r;
        end
        
        function test_TimeDependentOutput
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin);
            s.x0 = ts.x0;
            s.C = dscomponents.PointerOutputConv(@(t,mu)t,true);
            m.simulate();
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.simulate();
            clear s m r;
        end
    end
end
