classdef LocalLipKernelEstimator < error.BaseEstimator
    %COMPWISEKERNELESTIMATOR Summary of this class goes here
    %
    %   Requirements:
    %   ReducedModel.System.f is an instance of AKernelCoreFun
    %   FullModel.Approx is an instance of CompwiseKernelCoreFun
    %
    
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
        M1;
        M2;
        M3;
        Ma_norms;
        xi;
    end
    
    properties(Access=private,Transient)
        % Iteration stuff
        times;
        divals;
        e1vals;
        C;
        stepcnt;
        neg_e1=false;
    end
    
    methods
        function this = LocalLipKernelEstimator(rmodel)
            
            % Validity checks
            msg = error.LocalLipKernelEstimator.validModelForEstimator(rmodel);
            if ~isempty(msg)
                error(msg);
            end
            
            % Call superclass constructor with model argument
            this = this@error.BaseEstimator(rmodel);
            
            this.KernelLipschitzFcn = rmodel.System.f.SystemKernel.getLipschitzFunction;

            this.ExtraODEDims = 3;
            
            fm = rmodel.FullModel;
            % Obtain the correct snapshots
            if ~isempty(fm.Approx)
                % Standard case: the approx function is a kernel expansion.
                
                % Get centers of full approximation
                this.xi = fm.Approx.snData.xi;
                % Get full d x N coeff matrix of approx function
                Ma = fm.Approx.Ma;
            else
                % This is the also possible case that the full core
                % function of the system is a KernelExpansion.
                
                % Get centers of full core function
                this.xi = fm.System.f.snData.xi;
                % Get full d x N coeff matrix of core function
                Ma = fm.System.f.Ma;
            end
            if rmodel.System.f.RotationInvariantKernel
                this.xi = rmodel.System.f.snData.xi;
            end
            
            % Precompute norms
            this.Ma_norms = sqrt(sum(Ma.^2));
            
            % Only prepare matrices if projection was used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                % P = (I-VW^t)
                P = (eye(size(rmodel.V,1)) - rmodel.V*rmodel.W');
                D = P' * rmodel.G * P;
                
                hlp = Ma'*D*Ma;
                this.M1 = hlp;
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(rmodel.FullModel.System.B)
                    try
                        B = rmodel.FullModel.System.B.evaluate([],[]);
                    catch ME%#ok
                        B = rmodel.FullModel.System.B.evaluate(0,rmodel.System.getRandomParam);
                        warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                    end
                    % Filter too small entries
                    hlp = B'*D*B;
                    this.M2 = hlp;
                    hlp = Ma'*D*B;
                    this.M3 = hlp;
                end
            else
                % No projection means no projection error!
                this.M1 = 0;
                this.M2 = 0;
                this.M3 = 0;
            end
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            % extract current error
            eold = x(end-this.ExtraODEDims+1:end);
            e = zeros(3,1);
            
            % Compute \alpha(t)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-this.ExtraODEDims), t, mu);
            % An input function u is set
            if nargin == 5
                e(1) = phi'*this.M1*phi + ut'*this.M2*ut + phi'*this.M3*ut;
                % No input case
            else
                e(1) = phi'*this.M1*phi;
            end
            if ~this.neg_e1 && e(1) < 1
                this.neg_e1 = true;
            end
            %e(1) = sqrt(abs(e(1)));
            e(1) = sqrt(max(e(1),0));
            
            % Compute the local lipschitz constant estimations
            if this.ReducedModel.System.f.RotationInvariantKernel
                z = x(1:end-this.ExtraODEDims);
            else
                z = this.ReducedModel.V*x(1:end-this.ExtraODEDims);
            end
            di = this.xi - repmat(z,1,size(this.xi,2));
            di = sqrt(sum(di.^2,1));
            
            %% Verbose debug block
%             k = this.ReducedModel.System.f.SystemKernel;
%             ts = tic;
%             ci1 = k.getLocalGradientLipschitz(di,Inf,t,mu);
%             beta1 = this.Ma_norms * ci1';
%             t = toc(ts);
%             fprintf('LGL beta: %f, c-time: %f\n', beta1, t);
%             
%             ts = tic;
%             ci2 = k.getLocalSecantLipschitz(di,Inf,t,mu);
%             beta2 = this.Ma_norms * ci2';
%             t = toc(ts);
%             fprintf('LSL beta: %f, c-time: %f\n', beta2, t);
%             
%             ts = tic;
%             ci3 = k.getImprovedLocalSecantLipschitz(di,Inf,t,mu);
%             beta3 = this.Ma_norms * ci3';
%             t = toc(ts);
%             fprintf('ILSL beta: %f, c-time: %f\n', beta3, t);
%             beta = min([beta1,beta2,beta3]);

            %% Normal computations
            
            % Standard (worst-) Case
            Ct = Inf;
            % Time-discrete computation
            if this.UseTimeDiscreteC
                Ct = eold(1) + exp(eold(2)).*eold(3);
            end
            beta = this.Ma_norms * this.KernelLipschitzFcn(di,Ct,t,mu)';
            e(2) = beta;
            e(3) = eold(1)*beta*exp(-eold(2));
            
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
            this.LastError = eint(1,:) + exp(eint(2,:)).*eint(3,:);
            
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
        end
        
        function e0 = getE0(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = [this.ReducedModel.getExo(mu); 0; 0];
        end 
        
        function clear(this)
            clear@error.BaseEstimator(this);
            this.times = [];
            this.e1vals = [];
            this.divals = [];
        end
        
        %% Getter & Setter
        function set.Iterations(this, value)
            if value > 0 && isa(this.ReducedModel.ODESolver,'solvers.MLWrapper')%#ok
                warning('errorEst:LocalLipEst',...
                    'Build-In Matlab solvers cannot be use with this Error Estimator if Iterations are turned on.\nSetting Iterations = 0.');
                this.Iterations = 0;
            end
            this.Iterations = value;
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
            this.C = this.LastError;
            for it = 1:this.Iterations
%                 fprintf('%s,%s, Iteration %d\n',cl.Name,func2str(this.KernelLipschitzFcn),it);
                
                % Set step counter to one
                this.stepcnt = 1;
                % Expand the C error to fit the times
                this.C = this.C(idx);
                [ti, eint] = solver.solve(odefun, t, this.getE0(mu));
                this.C = eint(1,:) + exp(eint(2,:)).*eint(3,:);         
            end
            this.LastError = this.C;
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
            
            e = zeros(3,1);
            e(1) = this.e1vals(idx);
            beta = this.Ma_norms * this.KernelLipschitzFcn(...
                this.divals(idx,:),this.C(idx),t,mu)';
            e(2) = beta;
            e(3) = eold(1)*beta*exp(-eold(2));
            
            this.stepcnt = this.stepcnt+1;
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = [];
            if ~isa(rmodel.System.f,'dscomponents.AKernelCoreFun')
                errmsg = 'The reduced model''s core function must be a subclass of dscomponents.AKernelCoreFun for this error estimator.'; 
            elseif ~isempty(rmodel.FullModel.Approx) 
                if ~isa(rmodel.FullModel.Approx,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'The full model''s approx function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
                end
            elseif ~isa(rmodel.FullModel.System.f,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'If no approximation is used, the full model''s core function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
            end
        end
        
        function res = test_LocalLipKernelEstimator
            res = true;
            m = models.synth.KernelTest(10);
            m.offlineGenerations;
            r = m.buildReducedModel;
            r.ErrorEstimator = error.LocalLipKernelEstimator(r);
            
%             try
%                 m.ODESolver = solvers.MLWrapper(@ode23);
%                 r.ErrorEstimator = error.LocalLipKernelEstimator(r);
%                 r.ErrorEstimator.Iterations = 1;
%             catch ME%#ok
%                 res = true;
%             end
            
            m.ODESolver = solvers.Heun;
            r.ErrorEstimator = error.LocalLipKernelEstimator(r);
            r.ErrorEstimator.Iterations = 4;
            
            [t,y] = r.simulate;%#ok
        end
    end
end
