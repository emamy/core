classdef BaseFullModel < models.BaseModel & IParallelizable
    %BASEFULLMODEL The base class for any KerMor detailed model
    %   Implementers of custom models are to inherit from this base class
    %   in order for it to work with KerMor.
    %   For custom models, the properties of this class (combined with
    %   those from BaseModel) can be set to influence the model behaviour
    %   and reduction methods.
    %   For the implementation of custom dynamical systems, refer to
    %   BaseDynSystem.
    %
    % @todo build in time-tracking for offline phase etc
    %
    % See also: models BaseModel BaseDynSystem
    %
    % @author Daniel Wirtz @date 16.03.2010
    
    properties
        % The sampling strategy the Model uses
        %
        % See also: sampling BaseSampler
        Sampler;
        
        % The full model's data container.
        % Defaults to an empty container.
        %
        % See also: ModelData
        Data;
        
        % The general.PODFixspace instance used to compute the new snapshot
        % vectors for a given trajectory.
        %
        % The default is to use the mode 'abs' and value '1' for a new
        % snaphsot from a given trajectory.
        %
        % See also: general.POD general.Orthonormalizer
        %PODFix;
        
        % The reduction algorithm for subspaces
        %
        % See also: spacereduction BaseSpaceReducer
        SpaceReducer;
        
        % The approximation method for the CoreFunction
        %
        % Defaults to scalar kernel SVR regression
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
        preApproximationTrainingCallback;
        
        % Advanced property. 
        % Must be a function handle taking the current model instance.
        %
        % The handle is invoked at the end of the off4_computeApproximation
        % method.
        %
        % See preApproximationTrainingCallback for details.
        postApproximationTrainingCallback;
    end
    
%     properties(Access=private)
%         phase = [0 0 0 0 0];
%     end
    
    methods
        
        function this = BaseFullModel
            % Creates a new instance of a full model.
            
            % Setting default values for properties that are handle classes
            % will have the negative side-effect of having each instance of
            % BaseFullModel initialized with the SAME instance of the
            % Sampler, Approx etc. instances which of course is NOT
            % desireable.
            % See
            % http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_oop/bsdtcz7.html#bsdu1g9-1
            % for details.
            
            % Setup default values for the full model's components
            this.Sampler = sampling.RandomSampler;
            this.SpaceReducer = spacereduction.PODReducer;
            %             p = general.PODFixspace;
            %             p.Mode = 'abs';
            %             p.Value = 1;
            %             this.PODFix = p;
            this.Approx = approx.CompWiseInt;
            this.Data = models.ModelData;
            this.System = models.BaseDynSystem;
        end
        
        function time = off1_createParamSamples(this)
            % Offline phase 1: Sample generation.
            
            % Sampling
            time = tic;
            if ~isempty(this.Sampler)
                this.Data.ParamSamples = this.Sampler.generateSamples(this);
            end
            time = toc(time);
        end
        
        function time = off2_genProjectionTrainData(this)
            % Offline phase 2: Snapshot generation for subspace computation
            %
            % @todo Optimize snapshot array augmentation by preallocation
            % (later will be some storage class)
            time = tic;
            
            num_s = max(1,this.Data.SampleCount);
            num_in = max(1,this.System.InputCount);
            
            % Compute system dimension using x0.
            mu = this.Data.getParams(1);
            dims = length(this.System.x0(mu));
            
            trajlen = length(this.Times);
            
            if num_s*num_in < matlabpool('size')
                % @TODO: switch back if changed!
                this.ComputeParallel = false;
            end
            
            %% Parallel - computation
            if this.ComputeParallel
                idxmat = general.Utils.createCombinations(1:num_s,1:num_in);
                sn = zeros(dims+3,trajlen,size(idxmat,2));
                
                fprintf('Starting parallel projection training data computation of %d trajectories on %d workers...\n',size(idxmat,2),matlabpool('size'));
                % Iterate through all param/input combinations
                paridx = idxmat(1,:); % Use sliced variables here
                inidx = idxmat(2,:);
                clear idxmat;
                
                parfor idx = 1:size(paridx,2)
                    
                    % Check for parameters
                    mu = []; munum = 0;
                    if this.Data.SampleCount > 0%#ok
                        munum = paridx(idx);
                        mu = this.Data.getParams(munum);
                    end
                    % Check for inputs
                    inputidx = []; innum = 0;
                    if this.System.InputCount > 0
                        inputidx = inidx(idx);
                        innum = inputidx;
                    end
                    
                    % Get trajectory
                    [t, x] = this.computeTrajectory(mu, inputidx);
                    
                    % Assign snapshot values
                    sn(:,:,idx) = [ones(1,trajlen)*munum; ones(1,trajlen)*innum; t; x];
                end
                this.Data.ProjTrainData = sn(:,:);
                
                %% Non-parallel computation
            else
                sn = zeros(dims+3,trajlen*num_s*num_in);
                try
                    wh = waitbar(0,'Initializing...');
                    
                    % Assume no parameters or inputs
                    mu = [];
                    munum = 0;
                    inputidx = [];
                    innum = 0;
                    
                    cnt = 0;
                    % Iterate through all input functions
                    for inidx = 1:num_in
                        % Iterate through all parameter samples
                        for pidx = 1:num_s
                            
                            % Display
                            perc = cnt/(num_in*num_s);
                            waitbar(perc,wh,sprintf('Generating projection training data... %2.0f%%',round(perc*100)));
                            
                            % Check for parameters
                            if this.Data.SampleCount > 0
                                mu = this.Data.getParams(pidx);
                                munum = pidx;
                            end
                            % Check for inputs
                            if this.System.InputCount > 0
                                inputidx = inidx;
                                innum = inidx;
                            end
                            
                            % Get trajectory
                            [t, x] = this.computeTrajectory(mu, inputidx);
                            
                            % Assign snapshot values
                            curpos = cnt*trajlen+1;
                            sn(:,curpos:curpos+trajlen-1) = ...
                                [ones(1,trajlen)*munum; ones(1,trajlen)*innum; t; x];
                            
                            cnt=cnt+1;
                        end
                    end
                    close(wh);
                catch ME
                    close(wh);
                    rethrow(ME);
                end
                this.Data.ProjTrainData = sn;
            end
            
            
            time = toc(time);
        end
        
        function time = off3_computeReducedSpace(this)
            % Offline phase 3: Generate state space reduction
            
            time = tic;
            
            % Clear before running, so that in case of errors the matrix
            % from old reductions is unset.
            this.Data.V = [];
            this.Data.W = [];
            if ~isempty(this.SpaceReducer)
                if size(this.Data.ProjTrainData,2) == 1
                    % Easy case: Source dimension is already one. Just set V = Id.
                    this.Data.V = 1;
                    this.Data.W = 1;
                    warning('KerMor:spacereduction',['System''s state dimension'...
                        'is already one; no effective reduction.']);
                end
                fprintf('Computing reduced space...\n');
                [this.Data.V this.Data.W] = this.SpaceReducer.generateReducedSpace(this);
            end
            
            time = toc(time);
        end
        
        function time = off4_genApproximationTrainData(this)
            % Generates the training data for the `\hat{f}`-approximation.
            %
            % @todo create IMatrixArgument-interface or so, to detect
            % whether one can evaluate the system function passing a matrix
            % argument instead of a vector. -> speedup
            time = tic;
            if ~isempty(this.Approx)
                
                % Select subset of projection training data
                atd = this.Approx.selectTrainingData(this.Data);
                
                % If projection is used, train approximating function in
                % centers projected into the subspace.
                if ~isempty(this.SpaceReducer)
                    atd(4:end,:) = (this.Data.V*this.Data.W')*atd(4:end,:);
                end
                this.Data.ApproxTrainData = atd;
                
                % Compute f-Values at training data
                % Deactivated due to immense overhead and matlab crashes.
                % investigate further
                if false && this.ComputeParallel
                    fval = zeros(size(this.Data.ProjTrainData,1)-3,size(atd,2));
                    atddata = atd(4:end,:);
                    pardata = atd(1,:);
                    if KerMor.App.Verbose > 0
                        fprintf('Starting parallel f-values computation at %d points on %d workers...\n',size(atd,2),matlabpool('size'));
                    end
                    parfor sidx=1:size(atd,2)
                        fval(:,sidx) = ...
                            this.System.f.evaluate(atddata(:,sidx),... % x
                            atd(3,sidx),... % t
                            this.Data.getParams(pardata(sidx))); %#ok mu
                    end
                    this.Data.ApproxfValues = fval;
                else
                    this.Data.ApproxfValues = zeros(size(this.Data.ProjTrainData,1)-3,size(atd,2));
                    if KerMor.App.Verbose > 0
                        fprintf('Serial computation of f-values at %d points ...\n',size(atd,2));
                    end
                    for sidx=1:size(atd,2)
                        this.Data.ApproxfValues(:,sidx) = ...
                            this.System.f.evaluate(atd(4:end,sidx),... % x
                            atd(3,sidx),... % t
                            this.Data.getParams(atd(1,sidx))); % mu
                    end
                end
                
            end
            time = toc(time);
        end
        
        function time = off5_computeApproximation(this)
            % Offline phase 5: Core function approximation
            %
            % @todo proper subset selection algorithm/method
            time = tic;
            
            if ~isempty(this.Approx)
                if isempty(this.Data.ApproxTrainData) || isempty(this.Data.ApproxfValues)
                    error('No approximation training data available. Called off4_genApproximationTrainData?');
                end
                
                if isa(this.preApproximationTrainingCallback,'function_handle')
                    this.preApproximationTrainingCallback();
                end
                
                if KerMor.App.Verbose > 0
                    fprintf('Serial approximation computation for %d dimensions ...\n',size(this.Data.ApproxfValues,1));
                end
                
                % Approximate core function (is parallelizable for its own)
                this.Approx.approximateCoreFun(this);
                
                if isa(this.postApproximationTrainingCallback,'function_handle')
                    this.postApproximationTrainingCallback();
                end
            end
            
            time = toc(time);
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
            times(2) = this.off2_genProjectionTrainData;
            times(3) = this.off3_computeReducedSpace;
            times(4) = this.off4_genApproximationTrainData;
            times(5) = this.off5_computeApproximation;
            
            % Set time dirt flag to false as current sn fit the
            % times used.
            this.TimeDirty = false;
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
            
%             if this.TimeDirty
%                 warning(['The T or dt parameters have been changed since the last offline generations.\n'...
%                     'A call to offlineGenerations is required.']);
%             else
            if isempty(this.Data) || isempty(this.Data.ProjTrainData)
                error('No Snapshot data available. Forgot to call offlineGenerations before?');
            end
            tic;
            reduced = models.ReducedModel(this);
            time = toc;
        end
        
        % Not longer used as complete trajectories aren't stored in the
        % Data.Snapshot array
        %         function [t,x] = computeTrajectory(this, mu, inputidx)
        %             % Overrides the base method in BaseModel. For speed reasons any
        %             % already computed trajectories are returned without being
        %             %
        %             % Parameters:
        %             % mu: The parameter `\mu` to use. Set [] for none.
        %             % inputidx: The input function `u(t)` index to use. Set [] for
        %             % none.
        %             %
        %             % See also: models.BaseModel#computeTrajectory(mu, inputidx);
        %             %
        %             % @todo Possibly save the later computed trajectories in the
        %             % ModelData class?
        %
        %             if ~isempty(this.Data) && ~isempty(this.Data.ProjTrainData)
        %                 x = this.Data.getTrajectory(mu,inputidx);
        %                 if ~isempty(x)
        %                     t = this.Times;
        %                     return;
        %                 end
        %             end
        %             [t,x] = computeTrajectory@models.BaseModel(this, mu, inputidx);
        %
        %             % HERE: Save new trajectory in ModelData!
        %         end
        
    end
    
    %     methods(Access=protected,Sealed)
    %         function opts = trajectoryCompInit(this, mu, inputidx)%#ok
    %             % Implements the template method from models.BaseModel
    %             %
    %             % Here in the full model nothing is to do yet as only error
    %             % computation for reduced systems is performed so far in
    %             % ReducedModel.
    %             %
    %             % See also: models.BaseModel models.ReducedModel
    %
    %             % Nothing to do here yet.
    %             opts = [];
    %         end
    %     end
    
    %% Getter & Setter
    methods
        function set.Sampler(this, value)
            this.checkType(value, 'sampling.BaseSampler');%#ok
            this.Sampler = value;
        end
        
        function set.Data(this, value)
            this.checkType(value, 'models.ModelData');%#ok
            this.Data = value;
        end
        
        function set.SpaceReducer(this, value)
            this.checkType(value, 'spacereduction.BaseSpaceReducer');%#ok
            this.SpaceReducer = value;
        end
        
        function set.Approx(this, value)
            this.checkType(value, 'approx.BaseApprox');%#ok
            this.Approx = value;
        end
    end
    
    methods(Static)
        function test_BaseModels
            m = models.BaseFullModel;
            af = dscomponents.AffLinCoreFun;
            af.addSummand(@(t,mu)1, rand(4,4));
            af.addSummand(@(t,mu)sum(mu)*5, rand(4,4)*3);
            m.System.f = af;
            m.System.x0 = @(mu)sin(1:4)';
            m.offlineGenerations;
            red = m.buildReducedModel;
            % Test simulations
            m.simulate();
            red.simulate();
            
            % dont forget to free resources! (static handles seem to
            % persist over several method calls)
            clear af m red;
        end
        
        function test_BareModel
            m = models.BaseFullModel;
            m.SpaceReducer = [];
            m.Approx = [];
            m.Sampler = [];
            A = rand(2,2);
            m.System.f = dscomponents.PointerCoreFun(@(x,t,mu)A*x);
            m.System.x0 = @(mu)sin(1:2);
            m.offlineGenerations;
            red = m.buildReducedModel;
            red.simulate();
            % dont forget to free resources! (static handles seem to
            % persist over several method calls)
            clear af m red;
        end
        
        function test_LinearModel
            ts = testing.testsettings;
            m = ts.m;
            s = m.System;
            s.f = dscomponents.PointerCoreFun(ts.flin);
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
            s.f = dscomponents.PointerCoreFun(ts.flin);
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
            s.f = dscomponents.PointerCoreFun(ts.flin);
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
            s.f = dscomponents.PointerCoreFun(ts.flin_p);
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
            s.f = dscomponents.PointerCoreFun(ts.flin);
            s.x0 = ts.x0;
            s.Inputs = ts.Inputs;
            m.simulate([],1);
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
            s.f = dscomponents.PointerCoreFun(ts.flin_p);
            s.B = dscomponents.PointerInputConv(ts.B_p);
            s.x0 = ts.x0_p;
            for idx = 1:length(ts.params)
                p = ts.params(idx);
                s.addParam(p.Name, [p.MinVal p.MaxVal], p.Desired);
            end
            s.Inputs = ts.Inputs;
            m.simulate(m.System.getRandomParam, 1);
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