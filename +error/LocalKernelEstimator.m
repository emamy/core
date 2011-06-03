classdef LocalKernelEstimator < error.BaseKernelEstimator
    % LocalKernelEstimator: A-posteriori error estimator for kernel-based systems using local
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
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
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
        % @type error.BaseLocalLipschitzFunction
        %
        % See also: kernels.BaseKernel kernels.BellFunction error.BaseLocalLipschitzFunction
        LocalLipschitzFcn;
    end
    
    properties(Dependent)
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
        % Defaults to true. (As is best estimator atm)
        UseTimeDiscreteC;
    end
    
    properties(Transient, SetAccess=private)
        % The `d_i(t)` values for each integration time-step `t`.
        d_iValues;
    end
    
    properties(Access=private)
        Ma_norms;
        xi;
        mu;
        lstPreSolve;
        lstPostSolve;
        fIterations = 0;
        fTDC = true;
    end
    
    properties(Transient, Access=private)
        % Iteration stuff
        errEst;
        tstep;
    end
    
    methods
        function this = LocalKernelEstimator(rmodel)
            this = this@error.BaseKernelEstimator;
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.LocalKernelEstimator;
            % ExtraODEDims is set in constructor!
            copy = clone@error.BaseKernelEstimator(this, copy);
            % Clone local lipschitz function
            copy.LocalLipschitzFcn = this.LocalLipschitzFcn.clone;
            copy.fIterations = this.fIterations;
            copy.fTDC = this.fTDC;
            copy.Ma_norms = this.Ma_norms;
            copy.xi = this.xi;
            copy.d_iValues = this.d_iValues;
            copy.errEst = this.errEst;
            copy.tstep = this.tstep;
            copy.mu = this.mu;
            copy.lstPreSolve = addlistener(copy.ReducedModel.ODESolver,'PreSolve',@copy.cbPreSolve);
            copy.lstPreSolve.Enabled = this.lstPreSolve.Enabled;
            copy.lstPostSolve = addlistener(copy.ReducedModel.ODESolver,'PostSolve',@copy.cbPostSolve);
            copy.lstPostSolve.Enabled = this.lstPostSolve.Enabled;
        end
        
        function setReducedModel(this, rmodel)
            % Overrides the setReducedModel method from error.BaseEstimator
            % and performs additional offline computations.
            
            % Call superclass method to perform standard estimator
            % computations
            setReducedModel@error.BaseKernelEstimator(this, rmodel);
            
            this.lstPreSolve = addlistener(this.ReducedModel.ODESolver,'PreSolve',@this.cbPreSolve);
            this.lstPreSolve.Enabled = false;
            this.lstPostSolve = addlistener(this.ReducedModel.ODESolver,'PostSolve',@this.cbPostSolve);
            this.lstPostSolve.Enabled = false;
            
            fm = this.ReducedModel.FullModel;
            % Obtain the correct snapshots
            if ~isempty(fm.Approx)
                % Standard case: the approx function is a kernel expansion.
                % Get centers of full approximation
                this.xi = fm.Approx.Centers.xi;
                % Precompute norms
                this.Ma_norms = fm.Approx.Ma_norms;
            else
                % This is the also possible case that the full core
                % function of the system is a KernelExpansion.
                %
                % Get centers of full core function
                this.xi = fm.System.f.Centers.xi;
                % Precompute norms
                this.Ma_norms = fm.System.f.Ma_norms;
            end
            if this.ReducedModel.System.f.RotationInvariant
                this.xi = this.ReducedModel.System.f.Centers.xi;
            end
            
            ker = this.ReducedModel.System.f.SystemKernel;
            % Assign Lipschitz function
            lfcn = error.ImprovedLocalSecantLipschitz(ker);
            % Pre-Compute bell function xfeats if applicable
            [x,X] = general.Utils.getBoundingBox(this.xi);
            d = norm(X-x);
            lfcn.precompMaxSecants(d*2,200);
            this.LocalLipschitzFcn = lfcn;
        end
        
        function clear(this)
            clear@error.BaseKernelEstimator(this);
            this.d_iValues = [];
        end
        
        function prepareConstants(this)
            % Call superclass method
            prepareConstants@error.BaseKernelEstimator(this);
            
            % Returns the initial error at `t=0` of the integral part.
            this.lstPreSolve.Enabled = true;
            this.lstPostSolve.Enabled = true;
            % Call ISimConstants update function to compute values that are constant during a
            % simulation.
            this.LocalLipschitzFcn.prepareConstants;
        end
        
    end
    
    %% Getter & Setter
    methods
        function set.LocalLipschitzFcn(this, value)
            if ~isa(value,'error.BaseLocalLipschitzFunction')
                error('LocalLipschitzFcn must be a error.BaseLocalLipschitzFunction subclass.');
            end
            this.LocalLipschitzFcn = value;
        end
        
        function set.Iterations(this, value)
            if ~isnonnegintscalar(value)
                error('Iterations value must be a non-negative integer.');
            end
            if this.fTDC && value > 0
                warning('error:LocalKernelEstimator:NoIterationsWithTD','Cannot use iterations in conjunction with time-discrete C(t). Value will have no effect.');
            end
            this.fIterations = value;
        end
        
        function set.UseTimeDiscreteC(this, value)
            if ~islogical(value)
                error('The value must be a logical');
            end
            if value && this.fIterations > 0
                warning('error:LocalKernelEstimator:NoIterationsWithTD','Cannot use iterations in conjunction with time-discrete C(t). Iteration value will have no effect now.');
            end
            this.fTDC = value;
        end
        
        function value = get.Iterations(this)
            value = this.fIterations;
        end
        
        function value = get.UseTimeDiscreteC(this)
            value = this.fTDC;
        end
    end
    
    methods(Access=protected)
        
        function b = getBeta(this, x, t, mu)
            % Compute the local lipschitz constant estimations
            if this.ReducedModel.System.f.RotationInvariant
                z = x(1:end-this.ExtraODEDims);
            else
                z = this.ReducedModel.V*x(1:end-this.ExtraODEDims);
            end
            di = this.xi - repmat(z,1,size(this.xi,2));
            di = sqrt(sum(di.^2,1));
            
            %% Normal computations
            % Standard (worst-) Case
            Ct = Inf;
            % Time-discrete computation
            if this.fTDC
                Ct = x(end);
                % Keep track of distances when iterations are used
            elseif this.fIterations > 0
                this.d_iValues(this.StepNr,:) = di;
            end
            b = this.Ma_norms * this.LocalLipschitzFcn.evaluate(di, Ct, t, mu)';
        end
        
        function postprocess(this, t, x, mu, inputidx)
            
            postprocess@error.BaseKernelEstimator(this, t, x, mu, inputidx);
            
            this.StateError = x(end,:);
            
            % PreSolve not needed during iterations
            this.lstPreSolve.Enabled = false;

            % Iteration stuff
            if ~this.fTDC && this.fIterations > 0
                % Switch on/off listeners
                
                this.mu = mu;
                solver = this.ReducedModel.ODESolver;
                e0 = this.getE0(mu);
                for it = 1:this.fIterations
                    % Set time-step counter to one
                    this.tstep = 1;
                    % Solve
                    [ti, e] = solver.solve(@this.iterationODEPart, t, e0);
                end
                this.StateError = e;
            end
            
            % PostSolve still needed during iterations
            this.lstPostSolve.Enabled = false;
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

%             try
            b = this.Ma_norms * this.LocalLipschitzFcn.evaluate(...
                this.d_iValues(idx,:), this.errEst(idx), t, this.mu)';
%             catch ME
%                 keyboard;
%             end
            e = b*eold + this.EstimationData(2,idx);
            
            this.EstimationData(3,this.tstep) = b;
            
            this.tstep = this.tstep+1;
        end
        
        function cbPreSolve(this, sender, data)%#ok
            %fprintf('cbPreSolve in local kernel est, Lipfun:%s, It:%d, TDC:%d\n',class(this.LocalLipschitzFcn),this.fIterations,this.fTDC);
            this.d_iValues = zeros(length(data.Times),size(this.xi,2));
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
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = validModelForEstimator@error.BaseKernelEstimator(rmodel);
            if isempty(errmsg) && ~isa(rmodel.System.f,'dscomponents.AKernelCoreFun')
                errmsg = 'The reduced model''s core function must be a subclass of dscomponents.AKernelCoreFun for this error estimator.'; 
            end
            if isempty(errmsg) && ~isa(rmodel.System.f.SystemKernel,'kernels.BellFunction')
                errmsg = 'The system''s kernel must be a kernels.BellFunction for this error estimator.';
            end
        end
        
        function res = test_LocalKernelEstimator
            res = true;
            m = models.synth.KernelTest(10);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.ErrorEstimator = error.LocalKernelEstimator(r);
            
%             try
%                 m.ODESolver = solvers.ode.sMLWrapper(@ode23);
%                 r.ErrorEstimator = error.LocalKernelEstimator(r);
%                 r.ErrorEstimator.Iterations = 1;
%             catch ME%#ok
%                 res = true;
%             end
            
            m.ODESolver = solvers.ode.Heun;
            r.ErrorEstimator = error.LocalKernelEstimator(r);
            r.ErrorEstimator.Iterations = 4;
            
            [t,y] = r.simulate;%#ok
        end
    end
    
%     methods(Static,Access=protected)
%         function s = loadobj(s)
%             s = loadobj@error.BaseKernelEstimator(s);
%             addlistener(s.ReducedModel.ODESolver,'PreSolve',@s.PreSolve);
%         end
%     end
end
