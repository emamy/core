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
    % @change(0,3,sa,2011-04-23) Implemented Setters for the properties LocalLipschitzFcn
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
        
        % Determines how many postprocessing iterations for the estimator
        % are performed.
        %
        % Has computationally no effect if UseTimeDiscreteC is switched on.
        %
        % Default: 0
        Iterations = 0;
                
        % For the local Lipschitz constant estimation the parameter C can
        % be chosen to equal the error from the last time step. This has to
        % be investigated more thoroughly as integration errors from the
        % solver may lead to a loss of rigorousity.
        %
        % Defaults to true. (As is best estimator atm)
        UseTimeDiscreteC = true;
    end
    
    properties(Access=private)
        Ma_norms;
        xi;
        mu;
    end
    
    properties(Access=private,Transient)
        % Iteration stuff
        divals;
        errEst;
        stepcnt;
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
            % Fcn does not need to be cloned.
            copy.LocalLipschitzFcn = this.LocalLipschitzFcn;
            copy.Iterations = this.Iterations;
            copy.UseTimeDiscreteC = this.UseTimeDiscreteC;
            copy.Ma_norms = this.Ma_norms;
            copy.xi = this.xi;
            copy.divals = this.divals;
            copy.errEst = this.errEst;
            copy.stepcnt = this.stepcnt;
        end
        
        function setReducedModel(this, rmodel)
            % Overrides the setReducedModel method from error.BaseEstimator
            % and performs additional offline computations.
            
            % Call superclass method to perform standard estimator
            % computations
            setReducedModel@error.BaseKernelEstimator(this, rmodel);
            
            fm = this.ReducedModel.FullModel;
            % Obtain the correct snapshots
            if ~isempty(fm.Approx)
                % Standard case: the approx function is a kernel expansion.
                % Get centers of full approximation
                this.xi = fm.Approx.Centers.xi;
                % Precompute norms
                this.Ma_norms = sqrt(sum(fm.Approx.Ma.^2,1));
            else
                % This is the also possible case that the full core
                % function of the system is a KernelExpansion.
                %
                % Get centers of full core function
                this.xi = fm.System.f.Centers.xi;
                % Precompute norms
                this.Ma_norms = sqrt(sum(fm.System.f.Ma.^2,1));
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
            this.divals = [];
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
            if value > 0 && (isa(this.ReducedModel.ODESolver,'solvers.ode.MLWrapper') || isa(this.ReducedModel.ODESolver,'solvers.ode.MLode15i'))%#ok
                warning('errorEst:LocalLipEst',...
                    'Build-In Matlab solvers cannot be use with this Error Estimator if Iterations are turned on.\nSetting Iterations = 0.');
                this.Iterations = 0;
            end
            this.Iterations = value;
        end
        
        function set.UseTimeDiscreteC(this, value)
            if ~islogical(value)
                error('The value must be a logical');
            end
            this.UseTimeDiscreteC = value;
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
            if this.UseTimeDiscreteC
                Ct = x(end);
            end
            b = this.Ma_norms * this.LocalLipschitzFcn.evaluate(di, Ct, t, mu)';
            
            % Iteration stuff
            if this.Iterations > 0
                this.divals(end+1,:) = di;
            end
        end
        
        function postprocess(this, t, x, mu, inputidx)%#ok
            this.StateError = x(end,:);

            % Iteration stuff
            if this.Iterations > 0
                %this.times(end+1) = t(end);
                if this.UseTimeDiscreteC
                    warning('error:LocalLipErrorEstimator','Using Iterations together with TimeDiscreteC will yield no improvement. Not performing iterations.');
                else
                    this.performIterations(t, mu);
                end
            end
        end
    end
    
    methods(Access=private)
        function performIterations(this, t, mu)
            solver = this.ReducedModel.ODESolver;
            
            % Validity checks
            ti = this.EstimationData(1,:);
            if any(abs(sort(ti)-ti) > 100*eps)
                error('Time values recorded are not monotoneously increasing. Cannot propertly extend C, aborting.');
            end
            
            % Computes an expansion index set idx.
            % Expands the currently given error at the possibly coarser
            % timesteps t towards the recorded timesteps times during
            % simulation. Makes a worst-case estimation by assuming the
            % greater error value for times between two successive t
            % values.
%             T = repmat(t,length(ti),1);
%             Times = repmat(ti',1,length(t));
%             idx = sum(T < Times,2)+1;
            
            this.mu = mu;
            this.errEst = this.StateError;
            for it = 1:this.Iterations                
                % Set step counter to one
                this.stepcnt = 1;
%                 % Expand the C error to fit the times
%                 this.errEst = this.errEst(idx);
                [ti, this.errEst] = solver.solve(@this.iterationODEPart, t, this.init(mu));
            end
            this.StateError = this.errEst;
        end
        
        function e = iterationODEPart(this, t, eold)
            % THIS command will only work if the same times are passed to
            % the ode function than in the previous computations. hence,
            % one cannot use blackbox-odesolvers like the matlab builtin
            % ones which may use different timesteps.
            %idx = find(this.times == t,1);
            idx = this.stepcnt;
            
            if abs(this.EstimationData(1,idx)-t) > 100*eps
                error('The ODE solver does not work as required by the iterative scheme.');
            end

            b = this.Ma_norms * this.LocalLipschitzFcn.evaluate(...
                this.divals(idx,:), this.errEst(idx), t, this.mu)';
            e = b*eold + this.EstimationData(2,idx);
            
            this.EstimationData(3,this.stepcnt) = b;
            
            this.stepcnt = this.stepcnt+1;
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
end
