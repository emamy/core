classdef IterationCompLemmaEstimator < error.BaseCompLemmaEstimator
    % IterationCompLemmaEstimator: A-posteriori error estimator for kernel-based systems using local
    % lipschitz constants.
    %
    % Implementation as in [WH10], but with updated ExtraODEDims and numerical computation.
    %
    % Model requirements:
    % ReducedModel.System.f is an instance of AKernelCoreFun
    % FullModel.Approx is an instance of CompwiseKernelCoreFun
    %
    % @author Daniel Wirtz @date 2010-08-10
    %
    % @change{0,5,dw,2011-07-07}
    % - New properties smallz_flag to correctly transform the input space
    % dimension on beta-computation depending on whether rotation invariant kernels are used or not.
    % - Also made the estimator work with kernels.ParamTimeKernelExpansions
    %
    % @change{0,5,dw,2011-07-04} Changed this class name to "IterationCompLemmaEstimator".
    %
    % @change{0,4,dw,2011-05-29}
    % - Changed this classes name to "LocalKernelEstimator".
    % - Restructured the error estimators to better adopt to the current
    % formulation. Now the KernelEstimators have a function getBeta instead of implementing the
    % evalODEPart by themselves.
    %
    % @change{0,4,dw,2011-05-25} Changed the computations to the comparison lemma type. This reduced
    % the needed extra ODE dimensions to one and speeds up the evaluation process.
    %
    % @change{0,4,dw,2011-05-23} Adopted to the new error.BaseEstimator interface with separate output
    % error computation.
    %
    % @change{0,3,dw,2011-05-02} Changed the implementation of the evalODEPart so that only two
    % extra ODE dimensions are needed. This avoids NaN entries when exponential values grow too big.
    %
    % @change{0,3,sa,2011-04-23} Implemented Setters for the properties LocalLipschitzFcn
    % and UseTimeDiscreteC
    %
    % @new{0,1,dw,2010-08-10} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties
        % The internal kernel Lipschitz function to use.
        %
        % Depending on the underlying kernel, a specific lipschitz constant
        % estimation can be used. Each kernel has the method
        % "getLipschitzFcn" which is by default implemented in BaseKernel.
        % The resulting function handle is used by default; however, to
        % enhance customizability, one can manually define a lipschitz
        % estimation function here. (For an example: see the BellFunction
        % kernel interface for different possibilities)
        %
        % Function handle signature:
        % di: The distances from the reduced system to the centers
        % C: The locality constant that determines how far away the true
        % system can maximally be (can be `\infty`)
        % t: The current time `t`
        % mu: The current parameter `\mu`
        %
        % @type error.lipfun.Base
        %
        % See also: kernels.BaseKernel kernels.BellFunction error.lipfun.Base
        LocalLipschitzFcn;
        
        % Determines how many postprocessing iterations for the estimator
        % are performed.
        %
        % Has computationally no effect if UseTimeDiscreteC is switched on.
        %
        % Default: 0
        Iterations;
        
        % For the local Lipschitz constant estimation the parameter C can
        % be chosen to equal the error from the last time step. This has to
        % be investigated more thoroughly as integration errors from the
        % solver may lead to a loss of rigorousity.
        %
        % @type logical @default true
        %
        % Defaults to true. (As is best estimator atm)
        UseTimeDiscreteC = true;
    end
    
    properties(Transient, SetAccess=private)
        % The `d_i(t)` values for each integration time-step `t`.
        d_iValues;
    end
    
    properties(Access=private)
        Ma_norms;
        c; % the expansion centers
        fIterations = 0;
        fTDC = true;
        G;
        smallz_flag = false;
    end
    
    properties(Transient, Access=private)
        % Iteration stuff
        errEst;
        tstep;
        lstPreSolve;
        lstPostSolve;
    end
    
    methods
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.IterationCompLemmaEstimator;
            % ExtraODEDims is set in constructor!
            copy = clone@error.BaseCompLemmaEstimator(this, copy);
            % Clone local lipschitz function
            copy.LocalLipschitzFcn = this.LocalLipschitzFcn.clone;
            copy.Iterations = this.Iterations;
            copy.UseTimeDiscreteC = this.UseTimeDiscreteC;
            copy.Ma_norms = this.Ma_norms;
            copy.c = this.c;
            copy.d_iValues = this.d_iValues;
            copy.errEst = this.errEst;
            copy.tstep = this.tstep;
            if ~isempty(copy.ReducedModel)
                copy.lstPreSolve = addlistener(copy.ReducedModel.ODESolver,'PreSolve',@copy.cbPreSolve);
                copy.lstPreSolve.Enabled = this.lstPreSolve.Enabled;
                copy.lstPostSolve = addlistener(copy.ReducedModel.ODESolver,'PostSolve',@copy.cbPostSolve);
                copy.lstPostSolve.Enabled = this.lstPostSolve.Enabled;
            end
            copy.G = this.G;
            copy.smallz_flag = this.smallz_flag;
        end
        
        function offlineComputations(this, fm)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            %
            % Parameters:
            % fm: The full model. @type models.BaseFullModel
            
            offlineComputations@error.BaseCompLemmaEstimator(this, fm);
            
            % Obtain the correct snapshots
            f = [];
            if ~isempty(fm.Approx)
                f = fm.Approx.Expansion;
            end
            if isempty(f)
                % This is the also possible case that the full core
                % function of the system is a KernelExpansion.
                %
                % Get centers of full core function
                f = fm.System.f.Expansion;
            end
            this.c = f.Centers;
            this.Ma_norms = f.Ma_norms;
            
            % Assign Lipschitz function
            lfcn = error.lipfun.ImprovedLocalSecantLipschitz(f.Kernel);
            % Pre-Compute bell function xfeats if applicable
            [x,X] = Utils.getBoundingBox(this.c.xi);
            d = norm(X-x);
            lfcn.precompMaxSecants(d*2,200);
            this.LocalLipschitzFcn = lfcn;
        end
        
        function prepared = prepareForReducedModel(this, rm)
            % Prepares this estimator for use with a given reduced model.
            % Basically uses the projection matrices and some other reduced quantities for
            % local storage.
            %
            % Parameters:
            % rm: The reduced model @type models.ReducedModel
            %
            % Return values:
            % prepared: A clone of this estimator, with accordingly projected components
            prepared = prepareForReducedModel@error.BaseCompLemmaEstimator(this, rm);
            
            prepared.lstPreSolve = addlistener(rm.ODESolver,'PreSolve',@prepared.cbPreSolve);
            prepared.lstPostSolve = addlistener(rm.ODESolver,'PostSolve',@prepared.cbPostSolve);
            prepared.G = rm.G;
            
            f = rm.FullModel.Approx;
            if isempty(f)
                f = rm.FullModel.System.f;
            end
            if f.Expansion.Kernel.IsRBF && ~isempty(rm.V)
                % Use the projected centers z_i from the reduces system with x_i = Vz_i in this
                % case.
                hlp = f.project(rm.V, rm.W);
                prepared.c = hlp.Expansion.Centers;
                clear hlp;
                % The norm is then ||Vz - x_i||_G = ||z-z_i||_V'GV
                prepared.G = rm.V'*(prepared.G*rm.V);
                % Set flag for small z state vectors to true
                prepared.smallz_flag = true;
            end
        end
        
        function b = getBeta(this, xfull, t)
            % Compute the local lipschitz constant estimations
            %
            % Parameters:
            % xfull: The current state variable vector @type colvec
            % t: The current time `t`
            % mu: The current parameter `\mu`
            
            % Project x variable back to full space if centers are not within the space spanned by V
            % (indicated by smallz_flag)
            x = xfull(1:end-1);
            if ~this.smallz_flag && ~isempty(this.ReducedModel.V)
                x = this.ReducedModel.V*x;
            end
            di = this.c.xi - repmat(x,1,size(this.c.xi,2));
            di = sqrt(sum(di.*(this.G*di),1));
            
            %% Normal computations
            % Standard (worst-) Case
            Ct = Inf;
            % Time-discrete computation
            if this.UseTimeDiscreteC
                Ct = xfull(end);
                % Keep track of distances when iterations are used
            elseif this.Iterations > 0
                this.d_iValues(this.StepNr,:) = di;
            end
            
            % Check if time and param kernels are used
            if (~isempty(this.mu))
                f = this.ReducedModel.System.f;
                hlp = this.Ma_norms ...
                    .* f.Expansion.TimeKernel.evaluate(t, this.c.ti) ...
                    .* f.Expansion.ParamKernel.evaluate(this.mu, this.c.mui);
            else
                hlp = this.Ma_norms;
            end
            
            b = hlp * this.LocalLipschitzFcn.evaluate(di, Ct)';
        end
        
        function clear(this)
            clear@error.BaseCompLemmaEstimator(this);
            this.d_iValues = [];
        end
        
        function ct = prepareConstants(this, mu, inputidx)
            % Return values:
            % ct: The time needed for preprocessing @type double
            
            if this.Iterations > 0 && this.UseTimeDiscreteC
                warning('Ambiguous configuration. Having Iterations and UseTimeDiscreteC set; preferring UseTimeDiscreteC');
            end
            
            % Call superclass method
            ct = prepareConstants@error.BaseCompLemmaEstimator(this, mu, inputidx);
            st = tic;
            % Returns the initial error at `t=0` of the integral part.
            if isempty(this.lstPreSolve)
                this.lstPreSolve = addlistener(this.ReducedModel.ODESolver,'PreSolve',@this.cbPreSolve);
            end
            this.lstPreSolve.Enabled = true;
            if isempty(this.lstPostSolve)
                this.lstPostSolve = addlistener(this.ReducedModel.ODESolver,'PostSolve',@this.cbPostSolve);
            end
            this.lstPostSolve.Enabled = true;
            % Call ISimConstants update function to compute values that are constant during a
            % simulation.
            this.LocalLipschitzFcn.prepareConstants;
            ct = ct + toc(st);
        end
        
        function ct = postProcess(this, x, t, inputidx)
            % Return values:
            % ct: The time needed for postprocessing @type double
            st = tic;
            this.StateError = x(end,:);
            
            % PreSolve not needed during iterations
            this.lstPreSolve.Enabled = false;
            
            % Iteration stuff
            if ~this.UseTimeDiscreteC && this.Iterations > 0
                % Switch on/off listeners
                
                solver = this.ReducedModel.ODESolver;
                e0 = this.getE0(this.mu);
                for it = 1:this.Iterations
                    % Set time-step counter to one
                    this.tstep = 1;
                    % Solve
                    [~, e] = solver.solve(@this.iterationODEPart, t, e0);
                end
                this.StateError = e;
            end
            
            % PostSolve still needed during iterations
            this.lstPostSolve.Enabled = false;
            
            ct = postProcess@error.BaseCompLemmaEstimator(this, x, t, inputidx);
            
            ct = ct + toc(st);
        end
        
    end
    
    %% Getter & Setter
    methods
        function set.LocalLipschitzFcn(this, value)
            if ~isa(value,'error.lipfun.Base')
                error('LocalLipschitzFcn must be a error.lipfun.Base subclass.');
            end
            this.LocalLipschitzFcn = value;
        end
        
        function set.Iterations(this, value)
            if isempty(value)
                value = 0;
            elseif ~isscalar(value) || value < 0
                error('Iterations value must be a non-negative integer.');
            end
            this.Iterations = value;
        end
        
        function set.UseTimeDiscreteC(this, value)
            if isempty(value)
                value = false;
            elseif ~islogical(value)
                error('The value must be a logical');
            end
            this.UseTimeDiscreteC = value;
        end
    end
    
    methods(Access=private)
        
        function e = iterationODEPart(this, t, eold)
            % THIS command will only work if the same times are passed to
            % the ode function than in the previous computations. hence,
            % one cannot use blackbox-odesolvers like the matlab builtin
            % ones which may use different timesteps.
            %idx = find(this.times == t,1);
            idx = this.tstep;
            
            if abs(this.EstimationData(1,idx)-t) > 100*eps
                error('The ODE solver does not work as required by the iterative scheme.');
            end
            
            if (~isempty(this.mu))
                f = this.ReducedModel.System.f;
                hlp = this.Ma_norms ...
                    .* f.Expansion.TimeKernel.evaluate(t, this.c.ti) ...
                    .* f.Expansion.ParamKernel.evaluate(this.mu, this.c.mui);
            else
                hlp = this.Ma_norms;
            end
            b = hlp * this.LocalLipschitzFcn.evaluate(...
                this.d_iValues(idx,:), this.errEst(idx))';
            
            e = b*eold + this.EstimationData(2,idx);
            
            this.EstimationData(3,this.tstep) = b;
            
            this.tstep = this.tstep+1;
        end
        
        function cbPreSolve(this, sender, data)%#ok
            %fprintf('cbPreSolve in local kernel est, Lipfun:%s, It:%d, TDC:%d\n',class(this.LocalLipschitzFcn),this.Iterations,this.UseTimeDiscreteC);
            this.d_iValues = zeros(length(data.Times),size(this.c.xi,2));
        end
        
        function cbPostSolve(this, sender, data)%#ok
            % This is a bit dangerous as this gets called both during full trajectory simulations
            % and iterations simulations. For the first case, a full d x n matrix is in data.States
            % and for the latter it's just the error estimation ode. However, both will work at this
            % place as (end,:) grabs the correct values either way.
            this.errEst = data.States(end,:);
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(model)
            % Validations
            errmsg = validModelForEstimator@error.BaseCompLemmaEstimator(model);
            if isempty(errmsg) && ~isempty(model.Approx) && ~isa(model.Approx,'dscomponents.ParamTimeKernelCoreFun')
                errmsg = 'The model''s approximation function must be a subclass of dscomponents.ParamTimeKernelCoreFun for this error estimator.';
            end
            if isempty(errmsg) && isa(model.System.f,'dscomponents.ParamTimeKernelCoreFun') ...
                    && ~isa(model.System.f.Expansion.Kernel,'kernels.BellFunction')
                errmsg = 'The system''s kernel must be a kernels.BellFunction for this error estimator.';
            end
        end
        
        function res = test_IterationCompLemmaEstimator
            m = models.synth.KernelTest(10);
            e = error.IterationCompLemmaEstimator;
            e.UseTimeDiscreteC = false;
            e.Iterations = 4;
            
            m.ErrorEstimator = e;
            m.offlineGenerations;
            r = m.buildReducedModel;
            
            r.simulate(r.getRandomParam,[]);
            e.UseTimeDiscreteC = true;
            r.simulate(r.getRandomParam,[]);
            
            res = true;
        end
    end
    
    %     methods(Static,Access=protected)
    %         function s = loadobj(s)
    %             s = loadobj@error.BaseCompLemmaEstimator(s);
    %             this.lstPreSolve = addlistener(s.ReducedModel.ODESolver,'PreSolve',@s.cbPreSolve);
    %             this.lstPreSolve.Enabled = false;
    %             this.lstPostSolve = addlistener(s.ReducedModel.ODESolver,'PostSolve',@s.cbPostSolve);
    %             this.lstPostSolve.Enabled = false;
    %         end
    %     end
end
