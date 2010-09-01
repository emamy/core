classdef BaseFullModel < models.BaseModel
    %BASEFULLMODEL The base class for any KerMor detailed model
    %   Implementers of custom models are to inherit from this base class
    %   in order for it to work with KerMor.
    %   For custom models, the properties of this class (combined with
    %   those from BaseModel) can be set to influence the model behaviour
    %   and reduction methods.
    %   For the implementation of custom dynamical systems, refer to
    %   BaseDynSystem.
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
        PODFix;
        
        % The reduction algorithm for subspaces
        %
        % See also: spacereduction BaseSpaceReducer
        SpaceReducer;
        
        % The approximation method for the CoreFunction
        %
        % Defaults to scalar kernel SVR regression
        Approx;
    end
    
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
            p = general.PODFixspace;
            p.Mode = 'abs';
            p.Value = 1;
            this.PODFix = p;
            this.Approx = approx.CompWiseSVR;
            this.Data = models.ModelData;
            this.System = models.BaseDynSystem;
        end
        
        function off1_generateSamples(this)
            % Offline phase 1: Sample generation.
            
            % Sampling
            if ~isempty(this.Sampler)
                this.Data.ParamSamples = this.Sampler.generateSamples(this);
            end
        end
        
        function off2_generateSnapshots(this)
            % Offline phase 2: Snapshot generation.
            %
            % @todo Optimize snapshot array augmentation by preallocation
            % (later will be some storage class)
            
            num_samples = this.Data.SampleCount;
            num_inputs = this.System.InputCount;
            
            % Compute system dimension using x0.
            mu = [];
            if num_samples > 0
                mu = zeros(this.System.ParamCount,1);
            end
            dims = length(this.System.x0(mu));
            
            snapshots = zeros(dims+3,0);
            f_vals = zeros(dims,0);
            
            try
                wh = waitbar(0,'Initializing snapshot generation...');
                cnt = 0;
                fixvec = zeros(dims+1,0);
                
                % Assume no parameters or inputs
                mu = [];
                munum = 0;
                inputidx = [];
                innum = 0;
                
                % Iterate through all input functions
                for inidx = 1:max(1,num_inputs)
                    % Iterate through all parameter samples
                    for pidx = 1:max(1,num_samples)
                        
                        % Display
                        perc = cnt/((num_inputs+1)*(num_samples+1));
                        waitbar(perc,wh,sprintf('Generating snapshots ... %2.0f %%',perc*100));
                        
                        % Check for parameters
                        if num_samples > 0
                            mu = this.Data.ParamSamples(:,pidx);
                            munum = pidx;
                        end
                        % Check for inputs
                        if num_inputs > 0
                            inputidx = inidx;
                            innum = inidx;
                        end
                        
                        % Get trajectory
                        [t, x] = this.computeTrajectory(mu, inputidx);
                        
                        % Compute basis extension
                        newx = this.PODFix.computePODFixspace([t; x], fixvec);
                        
                        % Catch special case where no parameters or inputs
                        % are used
                        if num_samples + num_inputs <= 1
                            % If the PODFixspace returns just one
                            % dimension, approximation etc. does not make
                            % sense.
                            if size(newx,2) == 1
                                maxlen = min(size(x,1)+1,length(this.Times));
                                fprintf('Note: No parameters or inputs are used and PODFixspace generated just one snapshot.\n Changing PODFixspace and recomputing using Mode=''abs'' and Value=%d\n',maxlen);
                                this.PODFix.Mode = 'abs';
                                this.PODFix.Value = maxlen;
                                newx = this.PODFix.computePODFixspace([t; x], fixvec);
                            elseif size(newx,2) < .4*numel(this.Times)
                                warning('BaseFullModel:SmallBase','Snapshot base is very small, refine time steps or adjust PODFixspace.');
                            end
                        end
                        
                        fixvec = [fixvec newx];%#ok
                        
                        newlen = size(newx,2);
                        % Assign snapshot values
                        snapshots = [snapshots [ones(1,newlen)*munum;...
                                                ones(1,newlen)*innum;...
                                                newx]];%#ok
                        
                        %% Evaluate f at those points
                        for sidx=1:newlen
                            fx = this.System.f.evaluate(newx(2:end,sidx),newx(1,sidx),mu);
                            f_vals(:,end+1) = fx;%#ok
                        end
                        
                        cnt=cnt+1;
                    end
                end
                close(wh);
            catch ME
                close(wh);
                rethrow(ME);
            end
            
            this.Data.Snapshots = snapshots;
            this.Data.fValues = f_vals;
        end
        
        function off3_generateReducedSpace(this)
            % Offline phase 3: Generate state space reduction
            
            % Clear before running, so that in case of errors the matrix
            % from old reductions is unset.
            this.Data.V = [];
            this.Data.W = [];
            if ~isempty(this.SpaceReducer)
                if size(this.Data.Snapshots,2) == 1
                    % Easy case: Source dimension is already one. Just set V = Id.
                    this.Data.V = 1;
                    this.Data.W = 1;
                    warning('KerMor:spacereduction',['System''s state dimension'...
                        'is already one; no effective reduction.']);
                end
                try
                    wh = waitbar(.5,'Computing reduced space...');
                    [this.Data.V this.Data.W] = this.SpaceReducer.generateReducedSpace(this);
                    close(wh);
                catch ME
                    close(wh);
                    rethrow(ME);
                end
            end
        end
        
        function off4_generateApproximation(this)
            % Offline phase 4: Core function approximation
            if ~isempty(this.Approx)
                this.Approx.approximateCoreFun(this);
            end
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
            
            this.off1_generateSamples;
            this.off2_generateSnapshots;
            this.off3_generateReducedSpace;
            this.off4_generateApproximation;
            
            % Set time dirt flag to false as current snapshots fit the
            % times used.
            this.TimeDirty = false;
        end
        
        function reduced = buildReducedModel(this)
            % Builds a reduced model from a full model.
            %
            % Before calling this method ensure that offlineGenerations was
            % called at least once to provide the model's necessary data
            % for the reduction process.
            %
            % See also: offlineGenerations
            % @docupdate
            
            if this.TimeDirty 
                warning(['The T or dt parameters have been changed since the last offline generations.\n'...
                       'A call to offlineGenerations is required.']);
            elseif isempty(this.Data) || isempty(this.Data.Snapshots)
                error('No Snapshot data available. Forgot to call offlineGenerations before?');
            end
            reduced = models.ReducedModel(this);
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
%             if ~isempty(this.Data) && ~isempty(this.Data.Snapshots)
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

