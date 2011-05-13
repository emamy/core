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
    
    properties
        % The full model's data container.
        % Defaults to an empty container.
        %
        % @propclass{data}
        %
        % See also: ModelData
        Data;
    end
    
    properties(SetObservable)
        % The sampling strategy the Model uses
        %
        % @propclass{important}
        %
        % See also: sampling BaseSampler
        Sampler;
        
        % The reduction algorithm for subspaces
        %
        % @propclass{important}
        %
        % @default spacereduction.PODReducer
        % See also: spacereduction BaseSpaceReducer
        SpaceReducer;
        
        % The approximation method for the CoreFunction
        %
        % @propclass{important}
        %
        % @default approx.AdaptiveCompWiseKernelApprox
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
        % @propclass{optional}
        %
        % @todo create class events from this
        preApproximationTrainingCallback;
        
        % Advanced property. 
        % Must be a function handle taking the current model instance.
        %
        % The handle is invoked at the end of the off4_computeApproximation
        % method.
        %
        % @propclass{optional}
        %
        % See also: preApproximationTrainingCallback
        postApproximationTrainingCallback;
    end
    
    properties(SetObservable, Dependent)
        % The indices of inputs to use for training data generation.
        %
        % @propclass{optional} Set subindices to restrict offline generation phase to a subset of
        % inputs.
        %
        % @default 1:InputCount
        TrainingInputs;
    end
    
    properties(SetAccess=private, Dependent)
        % Gets the number of inputs used for training.
        TrainingInputCount;
    end
    
    properties(Access=private)
        msg;
        pstats;
        fTrainingInputs = [];
    end
           
    methods
        
        function this = BaseFullModel
            % Creates a new instance of a full model.
            
            this = this@models.BaseModel;
            
            % Setup default values for the full model's components
            this.Sampler = sampling.RandomSampler;
            this.SpaceReducer = spacereduction.PODReducer;
            %             p = general.PODFixspace;
            %             p.Mode = 'abs';
            %             p.Value = 1;
            %             this.PODFix = p;
            this.Approx = approx.AdaptiveCompWiseKernelApprox;
            this.Data = models.ModelData;
            
            % Register default properties
            this.registerProps('Sampler','SpaceReducer','Approx',...
                'preApproximationTrainingCallback','postApproximationTrainingCallback',...
                'Data','TrainingInputs');
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
        
        function time = off2_genTrainingData(this)
            % Offline phase 2: Snapshot generation for subspace computation
            %
            % @todo 
            % - Optimize snapshot array augmentation by preallocation
            % (later will be some storage class)
            % - Remove waitbar and connect to messaging system
            time = tic;
            
            num_s = max(1,this.Data.SampleCount);
            num_in = max(1,this.TrainingInputCount);
            
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
                idxmat = general.Utils.createCombinations(1:num_s,this.System.TrainingInputs);
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
                    if this.TrainingInputCount > 0
                        inputidx = inidx(idx);
                        innum = inputidx;
                    end
                    
                    % Get trajectory
                    [t, x] = this.computeTrajectory(mu, inputidx);                    
                    
                    % Assign snapshot values
                    sn(:,:,idx) = [ones(1,trajlen)*munum; ones(1,trajlen)*innum; t; x];
                end
                this.Data.TrainingData = sn(:,:);
                
                %% Non-parallel computation
            else
                sn = zeros(dims+3,trajlen*num_s*num_in);

                % Assume no parameters or inputs
                mu = [];
                munum = 0;
                inputidx = [];
                innum = 0;

                if KerMor.App.Verbose > 0
                    fprintf('Generating projection training data...\n');
                    p = 0;
                end
                cnt = 0;
                % Iterate through all input functions
                for inidx = 1:num_in
                    % Iterate through all parameter samples
                    for pidx = 1:num_s

                        % Display
                        if KerMor.App.Verbose > 0
                            perc = cnt/(num_in*num_s);
                            if perc > p
                                fprintf('%2.0f%%\n',round(perc*100));
                                p = ceil(perc*10)/10;
                            end
                        end
                        
                        % Check for parameters
                        if this.Data.SampleCount > 0
                            mu = this.Data.getParams(pidx);
                            munum = pidx;
                        end
                        % Check for inputs
                        if this.TrainingInputCount > 0
                            inputidx = this.TrainingInputs(inidx);
                            innum = inputidx;
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
                this.Data.TrainingData = sn;
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
                if size(this.Data.TrainingData,2) == 1
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
            % Generates the training data `x_i` for the
            % `\hat{f}`-approximation and precomputes `f(x_i)` values.
            %
            % @todo 
            % - include/check MultiArgumentEvaluations-possibility in parallel execution code
            % - Deactivated due to immense overhead and matlab crashes. investigate further
            time = tic;
            if ~isempty(this.Approx)
                
                % Select subset of projection training data
                atd = this.Approx.TrainDataSelector.selectTrainingData(this);
                
                % If projection is used, train approximating function in
                % centers projected into the subspace.
                if ~isempty(this.Data.V) && ~isempty(this.Data.W)
                    atd(4:end,:) = this.Data.V*(this.Data.W'*atd(4:end,:));
                end
                this.Data.ApproxTrainData = atd;
                
                % Compute f-Values at training data
                if this.ComputeParallel
                    fval = zeros(size(this.Data.TrainingData,1)-3,size(atd,2));
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
                    this.Data.ApproxfValues = zeros(size(this.Data.TrainingData,1)-3,size(atd,2));
                    if KerMor.App.Verbose > 0
                        fprintf('Serial computation of f-values at %d points ...\n',size(atd,2));
                    end
                    this.Data.ApproxfValues = ...
                            this.System.f.evaluate(atd(4:end,:),... % x
                            atd(3,:),... % t
                            this.Data.getParams(atd(1,:))); % mu
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
            times(2) = this.off2_genTrainingData;
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
            if isempty(this.Data) || isempty(this.Data.TrainingData)
                error('No Snapshot data available. Forgot to call offlineGenerations before?');
            end
            tic;
            reduced = models.ReducedModel(this);
            time = toc;
        end
        
        function [t,x] = computeTrajectory(this, mu, inputidx)
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
            
            this.checkProperties;
            
            [t,x] = computeTrajectory@models.BaseModel(this, mu, inputidx);
        end
        
        function printPropertyChangedSummary(this,levels)
            % Prints a summary about the properties of different levels which have not been changed
            % from their default setting.
            %
            % Call this method after any calls to models.BaseModel.simulate or
            % models.BaseModel.getTrajectory
            %
            % Parameters:
            % levels: [Optional] The property levels to print reports for. A list of admissible
            % values can be obtained by KerMorObject.getPropClasses. Default is to print ''all''
            % data.
            if nargin < 2
                levels = KerMorObject.getPropClasses;
            elseif ischar(levels)
                levels = {levels};
            end
            if ~isempty(this.pstats)
                c = this.pstats;
                total = sum(c,2);
                total(3) = total(3)/size(c,2);
                
                col = [1-total(3)/100 total(3)/100 0];
                cprintf(col,'Total unchanged properties: %d of %d (%2.2f%%%%)\n',total);
                for lidx = 1:length(levels)
                    col = [1-c(3,lidx)/100 c(3,lidx)/100 0];
                    lvidx = find(strcmp(levels{lidx},KerMorObject.getPropClasses),1);
                    cprintf(col, 'Unchanged ''%s'': %d of %d (%2.2f%%%%)\n',levels{lidx},c(:,lvidx));
                end
            end
        end
        
        function printPropertyChangedReport(this,levels)
            % Prints a detailed report about the properties at each level which have not been
            % changed from their default setting.
            %
            % Call this method after any calls to models.BaseModel.simulate or
            % models.BaseModel.getTrajectory
            %
            % Parameters:
            % levels: [Optional] The property levels to print reports for. A list of admissible
            % values can be obtained by KerMorObject.getPropClasses. Default is to print ''all''
            % data.
            
            if nargin < 2
                levels = KerMorObject.getPropClasses;
            elseif ischar(levels)
                levels = {levels};
            end
            this.printPropertyChangedSummary(levels);
            if ~isempty(this.msg)
                for idx = 1:length(levels)
                    m = this.msg.(levels{idx});
                    if ~isempty(m)
                        cprintf([.5 .5 0],'Messages for property class %s:\n',levels{idx});
                        fprintf('%s\n',m{:});
                        fprintf('\n');
                    end
                end
            end
        end
    end
        
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
        
        function set.TrainingInputs(this, value)
            if ~isposintmat(value)
                error('Value may only contain valid indices for the Inputs cell array.');
            elseif any(value > this.System.InputCount) || any(value < 1)
                error('Invalid indices for Inputs.');
            end
            this.fTrainingInputs = value;
        end
        
        function ti = get.TrainingInputs(this)
            ti = this.fTrainingInputs;
            if isempty(ti) && ~isempty(this.System)
                ti = 1:this.System.InputCount;
            end
        end
        
        function c = get.TrainingInputCount(this)
            c = length(this.TrainingInputs);
        end
    end
    
    methods(Access=private)        
         function checkProperties(this)
             % Checks all the model's properties recursively for unchanged default settings
             counts = struct;
             notchanged = struct;
             messages = struct;
             levels = KerMorObject.getPropClasses;
             for lidx=1:length(levels)
                 counts.(levels{lidx}) = 0;
                 messages.(levels{lidx}) = {};
             end
             notchanged = counts;
             
             %% Run recursive check
             recurCheck(this, general.collections.Dictionary, this.Name);
             
             % Store collected messages
             this.msg = messages;
             
             % Some total stats now
             c = [struct2array(notchanged); struct2array(counts)];
             c(3,:) = round(10000 * c(1,:) ./ c(2,:))/100;
             c(3,isnan(c(3,:))) = 100;
             this.pstats = c;
             
             if KerMor.App.Verbose > 1
                 total = sum(c,2);
                 total(3) = total(3)/size(c,2);

                 col = [total(3)/100 1-total(3)/100 0];
                 cprintf(col,'Total unchanged properties: %d of %d (%2.2f%%%%)\n',total);
             end
             
             % Issue warning if some critical properties are still unchanged
             if notchanged.critical > 0
                 link = 'critical properties';
                 if ~isempty(this.WorkspaceVariableName)
                    link = sprintf('<a href="matlab:%s.printPropertyChangedReport(''critical'')">critical properties</a>',this.WorkspaceVariableName);
                 end
                 fprintf(['SIMULATION RESULTS QUESTIONABLE: %d of %d ' link ' are still at their default value.\n'],...
                     notchanged.critical,counts.critical); %'Type <modelvarname>.printPropertyChangedReport(''critical'') for a detailed report.\n'],...
             end
             
             function recurCheck(obj, done, lvl)
                 mc = metaclass(obj);
                 done(obj.ID) = true;
                 
                 %% Link name preparations
                 objlink = editLink(mc.Name);
                 
                 %% Check the local properties
                 pc = obj.PropertiesChanged;
                 for pidx = 1:length(mc.Properties)
                     p = mc.Properties{pidx};
                     if strcmp(p.GetAccess,'public') && ~p.Constant && ~p.Transient && ~strcmp(p.SetAccess,'private')
                         key = [p.DefiningClass.Name '.' p.Name];
                         if pc.containsKey(key)
                             p = pc(key);
                             counts.(p.Level) = counts.(p.Level) + 1;
                             if ~p.Changed
                                 notchanged.(p.Level) = notchanged.(p.Level) + 1;
                                 if ~any(strcmp(p.Level,'data'))
                                     hlp = messages.(p.Level);
                                     %if strcmp(mc.Name,p.)
                                     hlp{end+1} = sprintf('%s is still unchanged!\nProperty brief: %s\nPropclass tag description:\n%s\n',...
                                         [lvl '[' objlink '] -> ' p.Name],p.Short,p.Text);%#ok
                                     messages.(p.Level) = hlp;
                                 end
                             end
                         elseif ~p.SetObservable && ~pc.containsKey([p.DefiningClass.Name '.' p.Name])
                             link2 = editLink(p.DefiningClass.Name);
                             fprintf('WARNING: Property %s of class %s is not <a href="matlab:docsearch SetObservable">SetObservable</a> but a candidate for a user-definable public property!\nFor more details see <a href="%s/propclasses.html">Property classes and levels</a>\n\n',p.Name,link2,KerMor.App.DocumentationLocation);
                         end
                         pobj = obj.(p.Name);
                         % Recursively register subobject's properties
                         if isa(pobj, 'KerMorObject') && isempty(done(pobj.ID))
                             recurCheck(pobj, done, [lvl '[' objlink '] -> ' p.Name]);
                         end
                     end
                 end
             end
             
             function l = editLink(classname)
                 dotpos = strfind(classname,'.');
                 if ~isempty(dotpos)
                     lname = classname(dotpos(end)+1:end);
                 else
                     lname = classname;
                 end
                 l = sprintf('<a href="matlab:edit %s">%s</a>',classname,lname);
             end
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
            %m.Approx = approx.DefaultCompWiseKernelApprox;
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
            m.System.f = dscomponents.PointerCoreFun(@(x,t,mu)A*x);
            m.System.x0 = @(mu)sin(1:2);
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
