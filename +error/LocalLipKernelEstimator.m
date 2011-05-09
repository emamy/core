classdef LocalLipKernelEstimator < error.BaseLipKernelEstimator
    % LocalLipKernelEstimator: A-posteriori error estimator for kernel-based systems using local
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
    % @new{0,1,dw,2010-08-10} Added this class.
    %
    % @change{0,3,dw,2011-05-02} CLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mLocalLipKernelEstimator.mhanged the implementation of the evalODEPart so that only two
    % extra ODE dimensions are needed. This avoids NaN entries when exponential values grow too big.
    %
    % @change(0,3,sa,2011-04-23) Implemented Setters for the properties KernelLipschitzFcn
    % and UseTimeDiscreteC
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
        % See also: kernels.BaseKernel kernels.BellFunction
        KernelLipschitzFcn;
        
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
    end
    
    properties(Access=private,Transient)
        % Iteration stuff
        times;
        divals;
        e1vals;
        errEst;
        stepcnt;
        neg_e1=false;
    end
    
    methods
        function this = LocalLipKernelEstimator(rmodel)
            this.ExtraODEDims = 2;
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
        
        function copy = clone(this)
            % Creates a deep copy of this estimator instance.
            copy = error.LocalLipKernelEstimator;
            % ExtraODEDims is set in constructor!
            copy = clone@error.BaseLipKernelEstimator(this, copy);
            copy.KernelLipschitzFcn = this.KernelLipschitzFcn;
            copy.Iterations = this.Iterations;
            copy.UseTimeDiscreteC = this.UseTimeDiscreteC;
            copy.Ma_norms = this.Ma_norms;
            copy.xi = this.xi;
            copy.times = this.times;
            copy.divals = this.divals;
            copy.e1vals = this.e1vals;
            copy.errEst = this.errEst;
            copy.stepcnt = this.stepcnt;
            copy.neg_e1 = this.neg_e1;
        end
        
        function setReducedModel(this, rmodel)
            % Overrides the setReducedModel method from error.BaseEstimator
            % and performs additional offline computations.
            
            % Call superclass method to perform standard estimator
            % computations
            setReducedModel@error.BaseLipKernelEstimator(this, rmodel);
            
            this.KernelLipschitzFcn = this.ReducedModel.System.f.SystemKernel.getLipschitzFunction;
            
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
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            % extract current error
            eold = x(end-this.ExtraODEDims+1:end);
            e = zeros(this.ExtraODEDims,1);
            
            % Compute \alpha(t)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-this.ExtraODEDims), t, mu);
            % An input function u is set
            if nargin == 5
                e(1) = phi*this.M1*phi' + phi*this.M2*ut + ut'*this.M3*ut;
                % No input case
            else
                e(1) = phi*this.M1*phi';
            end
            if ~this.neg_e1 && e(1) < 0
                this.neg_e1 = true;
            end
            e(1) = sqrt(abs(e(1)));
            %e(1) = sqrt(max(e(1),0));
            
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
                Ct = eold(1) + eold(2);
            end
            beta = this.Ma_norms * this.KernelLipschitzFcn(di,Ct,t,mu)';
            e(2) = beta*(eold(1) + eold(2));
            
            % Iteration stuff
            if this.Iterations > 0
                this.times(end+1) = t;
                this.e1vals(end+1) = e(1);
                this.divals(end+1,:) = di;
            end
        end
        
        function process(this, t, x, mu, inputidx)%#ok
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('CompWiseErrorEstimator:process','Integral part is all zero. Attention!');
            end
            this.LastError = eint(1,:) + eint(2,:);
            
            if this.neg_e1
                disp('LocalLipKernelEstimator: Negative alpha(t) norms occurred. Used zero instead.');
                this.neg_e1 = false;
            end
            
            % Iteration stuff
            if this.Iterations > 0
                this.times(end+1) = t(end);
                if this.UseTimeDiscreteC
                    warning('error:LocalLipErrorEstimator','Using Iterations together with TimeDiscreteC will yield no improvement. Not performing iterations.');
                else
                    this.performIterations(t, mu);
                end
            end
            
            % Tranform to output error estimation
            C = this.ReducedModel.FullModel.System.C.evaluate(0,[]);
            this.LastError = norm(C)*this.LastError;
        end
        
        function e0 = getE0(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = [this.ReducedModel.getExo(mu); 0];
        end 
        
        function clear(this)
            clear@error.BaseEstimator(this);
            this.times = [];
            this.e1vals = [];
            this.divals = [];
        end
        
        function set.KernelLipschitzFcn(this, value)
            this.KernelLipschitzFcn = value;
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
    
    methods(Access=private)
        function performIterations(this, t, mu)
            solver = this.ReducedModel.ODESolver;
            
            % Validity checks
            if any(abs(sort(this.times)-this.times) > 100*eps)
                error('Time values recorded are not monotoneously increasing. Cannot propertly extend C, aborting.');
            end
            
            % Computes an expansion index set idx.
            % Expands the currently given error at the possibly coarser
            % timesteps t towards the recorded timesteps times during
            % simulation. Makes a worst-case estimation by assuming the
            % greater error value for times between two successive t
            % values.
            T = repmat(t,length(this.times),1);
            Times = repmat(this.times',1,length(t));
            idx = sum(T < Times,2)+1;
            
%            cl = metaclass(this);
            odefun = @(t,x)this.iterationODEPart(x,t,mu);
            this.errEst = this.LastError;
            for it = 1:this.Iterations
%                 fprintf('%s,%s, Iteration %d\n',cl.Name,func2str(this.KernelLipschitzFcn),it);
                
                % Set step counter to one
                this.stepcnt = 1;
                % Expand the C error to fit the times
                this.errEst = this.errEst(idx);
                [ti, eint] = solver.solve(odefun, t, this.getE0(mu));
                this.errEst = eint(1,:) + eint(2,:);         
            end
            this.LastError = this.errEst;
        end
        
        function e = iterationODEPart(this, eold, t, mu)
            
            % THIS command will only work if the same times are passed to
            % the ode function than in the previous computations. hence,
            % one cannot use blackbox-odesolvers like the matlab builtin
            % ones which may use different timesteps.
            %idx = find(this.times == t,1);
            idx = this.stepcnt;
            
            if abs(this.times(idx)-t) > 100*eps
                error('The ODE solver does not work as required by the iterative scheme.');
            end
            
            e = zeros(this.ExtraODEDims,1);
            e(1) = this.e1vals(idx);
            beta = this.Ma_norms * this.KernelLipschitzFcn(...
                this.divals(idx,:),this.errEst(idx),t,mu)';
            e(2) = beta*(eold(1) + eold(2));
            
            this.stepcnt = this.stepcnt+1;
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = validModelForEstimator@error.BaseLipKernelEstimator(rmodel);
            if isempty(errmsg) && ~isa(rmodel.System.f,'dscomponents.AKernelCoreFun')
                errmsg = 'The reduced model''s core function must be a subclass of dscomponents.AKernelCoreFun for this error estimator.'; 
            end
        end
        
        function res = test_LocalLipKernelEstimator
            res = true;
            m = models.synth.KernelTest(10);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.ErrorEstimator = error.LocalLipKernelEstimator(r);
            
%             try
%                 m.ODESolver = solvers.ode.sMLWrapper(@ode23);
%                 r.ErrorEstimator = error.LocalLipKernelEstimator(r);
%                 r.ErrorEstimator.Iterations = 1;
%             catch ME%#ok
%                 res = true;
%             end
            
            m.ODESolver = solvers.ode.Heun;
            r.ErrorEstimator = error.LocalLipKernelEstimator(r);
            r.ErrorEstimator.Iterations = 4;
            
            [t,y] = r.simulate;%#ok
        end
    end
end
