classdef TPWLLocalLipEstimator < error.BaseEstimator
% TPWLLocalLipEstimator: Local-Lipschitz estimation based error estimator for reduced models
% obtained using the TPWL Approx class.
%
% Details about the math can be obtained upon contacting the author.
%
% @author Daniel Wirtz @date 2011-05-09
%
% @new{0,4,dw,2011-05-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % Determines how many postprocessing iterations for the estimator
        % are performed.
        %
        % Has computationally no effect if UseTimeDiscreteC is switched on.
        %
        % Default: 0
%         Iterations = 0;
        
        % For the local Lipschitz constant estimation the parameter C can
        % be chosen to equal the error from the last time step. This has to
        % be investigated more thoroughly as integration errors from the
        % solver may lead to a loss of rigorousity.
        %
        % Defaults to true. (As is best estimator atm)
        UseTimeDiscreteC = true;
    end
    
    properties(Access=private)
        KernelLipschitzFcn;
        A;
        Ab;
        b;
        AB;
        bB;
        B;
        AV;
        neg_e1 = false;
        Ainorms;
        xi;
    end
    
    methods
        function this = TPWLLocalLipEstimator
            this.ExtraODEDims = 2;
        end
        
        function copy = clone(this)
            error('not yet implemented');
        end
    
        function p = prepareForReducedModel(this, rmodel)
            % Overrides the setReducedModel method from error.BaseEstimator
            % and performs additional offline computations.
            
            % Call superclass method to perform standard estimator
            % computations
            p = prepareForReducedModel@error.BaseEstimator(this, rmodel);
                        
            p.KernelLipschitzFcn = rmodel.System.f.GaussWeight.getLipschitzFunction;
            
            % Get full model
            fm = rmodel.FullModel;
            if ~isempty(fm.System.B)
                try
                    fB = fm.System.B.evaluate([],[]);
                catch ME%#ok
                    fB = fm.System.B.evaluate(0,rmodel.getRandomParam);
                    warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                end
            end
            
            p.xi = fm.Approx.xi;
            
            Ai = fm.Approx.Ai;
            bi = fm.Approx.bi;
            G = rmodel.GScaled;
            n = length(Ai);
            % Perform any offline computations/preparations
            % Only prepare matrices if projection is used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                % Compute projected versions (tilde..)
                for i=1:n
                    Ai{i} = (Ai{i} - rmodel.V*(rmodel.W'*Ai{i}))*rmodel.V;
                end
                bi = bi - rmodel.V*(rmodel.W'*bi);
                
                p.A = cell(n,n);
                p.Ab = cell(n,n);
                p.b = zeros(n,n);
                [x,y] = meshgrid(1:n);
                all = [reshape(x,1,[]); reshape(y,1,[])];
                for k = 1:n*n
                    i = all(1,k);
                    j = all(2,k);
                    p.A{i,j} = Ai{i}'*(G*Ai{j});
                    p.Ab{i,j} = Ai{i}'*(G*bi(:,j));
                    %p.b(i,j) = bi(:,i)'*G*bi(:,j);
                end
                p.b = bi'*(G*bi);
                
                if ~isempty(fB)
                    p.AB = cell(1,n);
                    for i=1:n
                        p.AB{i} = Ai{i}'*(G*fB);
                    end
                    p.bB = bi'*(G*fB);
                    p.B = fB'*(G*fB);
                end
                
                % Compute AV_i`s for \beta(s) comp
                p.AV = cell(1,n);
                p.Ainorms = zeros(1,n);
                for i=1:n
                    p.AV{i} = Ai{i}'*(G*Ai{i});
                    p.Ainorms(i) = norm(Ai{i});
                end
            else
                error('\todo');
%                 % No projection means no projection error!
%                 n = size(Ma,2);
%                 
%                 if ~isempty(B)
%                     b = size(B,2);
%                     p.M2 = zeros(n,b);
%                     p.M3 = zeros(b,b);
%                 end
            end
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            % Evaluates the auxiliary ode part for the TPWL estimator.
            %
            % @attention: NOT WORKING PROPERLY YET.
            %
            % Parameters:
            % x: The full extended state variable vector. Extended means that
            % the last @ref ExtraODEDims rows contain the error estimators own
            % data. If not used, implementers must take care to ditch those
            % values if any function evaluations are performed within the
            % integral part. @type colvec
            % t: The current time `t` @type double
            % mu: The current parameter `\mu` @type colvec
            % ut: The value of the input function `u(t)` if given, [] else.
            % @type double
            %
            % Return values:
            % e: The auxiliary ode part value.
            % extract current error
            
            eold = x(end-this.ExtraODEDims+1:end);
            e = zeros(this.ExtraODEDims,1);
            wf = this.ReducedModel.System.f.GaussWeight;
            
            % Get actual z(s)
            z = x(1:end-this.ExtraODEDims);
            
            % Get distances
            di = this.xi - repmat(this.ReducedModel.V*z,1,size(this.xi,2));
            di = sqrt(sum(di.^2,1));
            di(di == 0) = eps;
            w = wf.evaluateScalar(di/min(di));
            w = w / sum(w);
            idx = 1:length(w);
            % Use same criterion as also used in normal approx evaluation
            idx = idx(w > this.ReducedModel.System.f.MinWeightValue);
            m = length(idx);
            [x,y] = meshgrid(idx);
            all = [reshape(x,1,[]); reshape(y,1,[])];
            e(1) = 0;
            % Run semi-quadratic loop as many weights hopefully are zero
            for k = 1:m*m
                i = all(1,k);
                j = all(2,k);
                e(1) = e(1) + w(i)*w(j)*(z'*this.A{i,j}*z + 2*z'*this.Ab{i,j} + this.b(i,j));
            end
            
            if nargin == 5 % An input function u is set
                for i=1:m
                    e(1) = e(1) + w(i)*(z'*this.AB{i}*ut + this.bB(i,:)*ut);
                end
                e(1) = e(1) + ut'*this.B*ut;
            end
            
            if ~this.neg_e1 && e(1) < 0
                this.neg_e1 = true;
            end
            %e(1) = sqrt(abs(e(1)));
            e(1) = sqrt(max(e(1),0));
            
            %% Normal computations
            % Standard (worst-) Case
            Ct = Inf;
            % Time-discrete computation
            if this.UseTimeDiscreteC
                Ct = eold(1) + eold(2);
            end
            
            % Part one of beta
            beta = wf.evaluateScalar(max(0,di-Ct)) * this.Ainorms';
            % Part two of beta
            av = cellfun(@(A)sqrt(z'*A*z),this.AV);
            beta = beta + av * this.KernelLipschitzFcn(di,Ct,t,mu)';
            
            e(2) = beta*(eold(1) + eold(2));
            
%             % Iteration stuff
%             if this.Iterations > 0
%                 this.times(end+1) = t;
%                 this.e1vals(end+1) = e(1);
%                 this.divals(end+1,:) = di;
%             end
        end
        
        function process(this, t, x, mu, inputidx)%#ok
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('CompWiseErrorEstimator:process','Integral part is all zero. Attention!');
            end
            this.StateError = eint(1,:) + eint(2,:);
            
            if this.neg_e1
                disp('TPWLLocalLipEstimator: Negative alpha(t) norms occurred. Used zero instead.');
                this.neg_e1 = false;
            end
            
%             % Iteration stuff
%             if this.Iterations > 0
%                 this.times(end+1) = t(end);
%                 if this.UseTimeDiscreteC
%                     warning('error:TPWLLocalLipEstimator','Using Iterations together with TimeDiscreteC will yield no improvement. Not performing iterations.');
%                 else
%                     this.performIterations(t, mu);
%                 end
%             end
            
            % Tranform to output error estimation
            C = this.ReducedModel.FullModel.System.C.evaluate(0,[]);
            this.StateError = norm(C)*this.StateError;
        end
        
        function e0 = init(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = [this.ReducedModel.getExo(mu); 0];
        end 
        
        function clear(this)
            clear@error.BaseEstimator(this);
            this.neg_e1 = false;
%             this.times = [];
%             this.e1vals = [];
%             this.divals = [];
        end
        
%         function set.Iterations(this, value)
%             if value > 0 && (isa(this.ReducedModel.ODESolver,'solvers.MLWrapper') || isa(this.ReducedModel.ODESolver,'solvers.MLode15i'))%#ok
%                 warning('errorEst:LocalLipEst',...
%                     'Build-In Matlab solvers cannot be use with this Error Estimator if Iterations are turned on.\nSetting Iterations = 0.');
%                 this.Iterations = 0;
%             end
%             this.Iterations = value;
%         end
    end
    
   
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = '';
            if ~isa(rmodel.FullModel.System.C,'dscomponents.LinearOutputConv')
                errmsg = 'Local Lipschitz estimators work only for constant linear output conversion.';
            end
            if isempty(errmsg) && ~isa(rmodel.FullModel.Approx,'approx.TPWLApprox')
                errmsg = 'The Approx class of the full model must be a subclass of approx.TPWLApprox for this error estimator.'; 
            end
        end
        
%         function res = test_TPWLLocalLipEstimator
%             res = true;
%             m = models.synth.KernelTest(10);
%             m.offlineGenerations;
%             r = m.buildReducedModel;
%             r.ErrorEstimator = error.IterationCompLemmaEstimator(r);
%             
% %             try
% %                 m.ODESolver = solvers.sMLWrapper(@ode23);
% %                 r.ErrorEstimator = error.IterationCompLemmaEstimator(r);
% %                 r.ErrorEstimator.Iterations = 1;
% %             catch ME%#ok
% %                 res = true;
% %             end
%             
%             m.ODESolver = solvers.Heun;
%             r.ErrorEstimator = error.IterationCompLemmaEstimator(r);
%             r.ErrorEstimator.Iterations = 4;
%             
%             [t,y] = r.simulate;%#ok
%         end
    end
    
end