classdef DEIMEstimator < error.BaseEstimator
    % DEIMEstimator: A-posteriori error estimation for DEIM reduced models.
    %
    % @author Daniel Wirtz @date 2012-05-10
    %
    % @change{0,6,dw,2012-06-11} 
    % - Added support for different norm-inducing matrices `G`.
    % - New property UseTrueDEIMErr to enable use of the actual DEIM
    % approximation error within the alpha term computation (Experimental
    % use)
    %
    % @new{0,6,dw,2012-06-08}
    % - Introduced the property JacMatDEIMMaxOrder to control the jacobian
    % DEIM approximation quality
    % - Made the JacMDEIM property private
    % - Introduced flags for more expensive (comparison-purpose)
    % estimators
    % - Removed the ErrorOrder setting as the DEIMEstimator directly depends
    % on the Order setting of the reduced model's DEIM approximation.
    % Setting the DEIM order "remotely" from the error estimator introduced
    % more trouble than necessary
    % - Not re-using the previously computed JacMDEIM approximation upon
    % new runs of offlineGenerations (mess)
    % - Giving out a warning if no sufficient reduction by the similarity
    % transformation is achieved and disabling it for less than 20% reduction. 
    % - Using an event listener for the DEIM OrderUpdated event to
    % re-compute local error estimation matrices
    %
    % @change{0,6,dw,2012-05-26} Updated the computation to the new structure
    % of systems (A + f components) and replaced the previous (errorneous) one
    %
    % @new{0,6,dw,2012-05-10} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties
        % The maximum size of the similarity transformation.
        %
        % Set to empty if not wanted. The BaseFullModel sets this value at
        % offline generations phase
        %
        % @type integer @default 50
        %
        % See also: JacSimTransSize
        JacSimTransMaxSize = 50;
        
        % Maximum order for the DEIM approximation of the state space
        % jacobian.
        %
        % @type integer @default 50
        JacMatDEIMMaxOrder = 50;
        
        % The selector that chooses the full model's trajectory points that
        % are used for the MatrixDEIM approximation of the jacobian and the
        % partial similarity transformation.
        %
        % @type data.selection.ASelector @default
        % data.selection.DefaultSelector
        TrainDataSelector;
        
        % "Expensive version": Using the full system's jacobian for
        % logarithmic norm computation.
        % Included for easy error estimator comparison.
        UseFullJacobian = false;
        
        % "Even More Expensive version": Precompute the full solution and
        % use the true logarithmic lipschitz constant of `f` in
        % computations.
        % Included for easy error estimator comparison.
        UseTrueLogLipConst = false;
        
        % "Even More Expensive version 2": Precompute the full solution and
        % use the true logarithmic lipschitz constant of the jacobian of `f` in
        % computations.
        % Included for easy error estimator comparison.
        UseJacobianLogLipConst = false;
        
        % "Expensive" version that uses the true approximation error
        % between the DEIM approximation and the full system function
        % instead of the `M'` variant.
        UseTrueDEIMErr = false;
    end
    
    properties%(SetAccess=private)
        % The complete similarity transformation matrix of size
        % `d \times JacSimTransMaxSize`
        %
        % (Debug/testing use)
        %
        % @type matrix<double> @default []
        QFull = [];
        
        % The singular values computed by the SimTransPOD in the offline
        % stage.
        %
        % (Debug/testing use)
        %
        % @type rowvec<double> @default []
        QSingVals = [];
        
        % The MatrixDEIM instance used to approximate the state space
        % jacobian matrices
        %
        % The default MaxOrder setting is 50.
        %
        % @type general.MatrixDEIM
        JacMDEIM;
        
        M3 = [];
        M4 = [];
        M5 = [];
        M6 = [];
        M7 = [];
        M8 = [];
        M9 = [];
        M10 = [];
        M11 = [];
        M12 = [];
        Ah; % the \hat{A} component
        Bh; % the \hat{B} component
        
        % The precomputed logarithmic norms of the system's A components
        Aln;
        
        % Stores zi min- and max scaling
        scale;
        
        % The local logarithmic norm estimation
        kexp;
    end
    
    properties(Dependent)
        % The size of the partial similarity transform to apply to the
        % matrix DEIM approximated jacobians.
        %
        % Upper bounded by JacSimTransMaxSize and defaults to the ceiled
        % half of the maximum size (set during offlineComputations)
        %
        % @type integer @default ceil(JacSimTransMaxSize/2)
        %
        % See also: JacSimTransMaxSize
        JacSimTransSize;
        
        % The actual order of the DEIM approximation for the system's
        % jacobian.
        %
        % Upper bounded by JacMatDEIMMaxOrder and defaults to the ceiled
        % half of the maximum size (set during offlineComputations)
        %
        % @type integer @default ceil(JacMatDEIMMaxOrder/2)
        %
        % See also: JacMatDEIMMaxOrder
        JacMatDEIMOrder;
    end
    
    properties(Access=private)
        jstSize = 0;
        deim;
    end
    
    properties(Access=private, Transient)
        fullsol;
        uolst = [];
        silent = false;
    end
    
    methods
        function this = DEIMEstimator
            this = this@error.BaseEstimator;
            this.ExtraODEDims = 1;
            this.TrainDataSelector = data.selection.DefaultSelector;
        end
        
        function offlineComputations(this, fm)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            
            if KerMor.App.Verbose > 0
                fprintf('error.DEIMEstimator: Starting offline computations...\n');
            end
            
            % Call superclass
            offlineComputations@error.BaseEstimator(this, fm);
            
            % Precompute log norm of A component if existing
            this.Aln = [];
            fA = fm.System.A;
            if ~isempty(fA)
                % Precompute logarithmic norm of A(t,\mu)
                a = general.AffParamMatrix;
                if isa(fA, 'dscomponents.LinearCoreFun')
                    a.addMatrix('1',general.Utils.logNorm(fA.A));
                else
                    for i=1:fA.N
                        a.addMatrix(['abs(' fA.funStr{i} ')'],...
                            general.Utils.logNorm(fA.getMatrix(i)));
                    end
                end
                this.Aln = a;
            end
            
            % Jacobian matrix DEIM approx
            this.compute_jacmdeim(fm);
            
            % Similarity transformation
            %this.compute_simtrans(fm);

            % Compute approximation for local logarithmic norm
            %             this.computeLogNormApprox;
        end
        
        function setReducedModel(this, rmodel)
            setReducedModel@error.BaseEstimator(this, rmodel);
            newd = rmodel.System.f;
            
            % Perform projection of JacMDeim instance
            if ~isempty(rmodel.V)
                this.JacMDEIM = this.JacMDEIM.project(rmodel.V, rmodel.W);
            end
            
            % Precompute all the projected quantities that are independent of the DEIM error
            % orders
            fm = rmodel.FullModel;
            md = fm.Data;
            fA = fm.System.A;
            G = fm.G;
            A = [];
            if ~isempty(fA)
                hlp = fA*md.V;
                A = hlp - md.V*(md.W'*hlp); % (I-VW^T)AV
                this.M6 = A'*(G*A);
                this.Ah = A;
            else
                this.M6 = []; this.M7 = []; this.M8 = []; this.M9 = [];
            end
            B = [];
            if ~isempty(fm.System.B)
                B = fm.System.B - md.V*(md.W'*fm.System.B);
                this.M12 = B'*(G*B);
                this.Bh = B;
            else
                this.M9 = []; this.M10 = []; this.M11 = []; this.M12 = [];
            end
            if ~isempty(A) && ~isempty(B)
                if isa(B,'dscomponents.LinearInputConv')
                    this.M9 = 2*A'*(G*B.B);
                else
                    % Both inherit from AffParamMatrix, so they can be
                    % multiplied :-)
                    this.M9 = 2*A'*(G*B);
                end
            else
                this.M9 = [];
            end
            
            if isempty(this.uolst) || (~isempty(this.deim) && ...
                    this.deim ~= newd)
                delete(this.uolst);
                this.uolst = addlistener(newd,'OrderUpdated',@this.handleOrderUpdate);
                this.deim = newd;
                this.updateErrMatrices;
            end
            this.deim = newd;
        end
        
        function ct = prepareConstants(this, mu, inputidx)
            if this.deim.Order(2) == 0 && ~this.UseTrueDEIMErr
                warning('KerMor:DEIMEstimator','No DEIM error order set. Disabling error estimator.');
                this.Enabled = false;
                ct = 0;
                return;
            end
            % True log lip const comparison estimator
            if this.UseTrueLogLipConst || this.UseJacobianLogLipConst
                [~, this.fullsol, ct] = this.ReducedModel.FullModel.computeTrajectory(mu,inputidx);
            else
                ct = 0;
            end
            this.kexp = [];
        end
        
        function a = getAlpha(this, x, t, mu, ut)
            x = x(1:end-1);
            % More expensive variant, using the true DEIM approximation
            % error instead of the estimated ones using M,M' technique
            if this.UseTrueDEIMErr
                rm = this.ReducedModel;
                fs = rm.FullModel.System;
                I = speye(size(rm.V,1));
                A = fs.A.evaluate(rm.V*x,t,mu);
                a = A - rm.V*(rm.W'*A);
                a = a + fs.f.evaluate(rm.V*x,t,mu) ...
                    - rm.V*rm.System.f.evaluate(x,t,mu);
                if ~isempty(fs.B)
                    a = a + (I-rm.V*rm.W')*fs.B.evaluate(t,mu)*ut;
                end
                a = Norm.LG(a,rm.FullModel.G);
            else
                f = this.ReducedModel.System.f.f;
                v1 = f.evaluateComponentSet(1, x, t, mu);
                v2 = f.evaluateComponentSet(2, x, t, mu);

                % NOTE: The factors 2 before some summands are added at offline
                % stage!!
                a = v1'*this.M3*v1 - v1'*this.M4*v2 + v2'*this.M5*v2;
                if ~isempty(this.M6) % An time/param affine A is set
                    a = a + x'*this.M6.evaluate(x,t,mu) + v1'*this.M7.evaluate(x,t,mu) ...
                        - v2'*this.M8.evaluate(x,t,mu);
                    if ~isempty(ut)
                        a = a + x'*this.M9.evaluate(ut,t,mu);
                    end
                end
                if ~isempty(ut) % An input function u is set
                    a = a + v1'*this.M10.evaluate(t,mu)*ut - v2'*this.M11.evaluate(t,mu)*ut ...
                        + ut'*this.M12.evaluate(t,mu)*ut;
                end
                a = sqrt(a);
%                 %% Validation: direct computation (expensive)
%                 V = this.ReducedModel.V;
%                 fs = this.ReducedModel.FullModel.System;
%                 I = speye(size(V,1));
%                 hlp = fs.A.evaluate(V*x,t,mu);
%                 a_1 = hlp - V*(V'*hlp);
%                 a_2 = this.ReducedModel.System.f.M1 * v1;
%                 a_3 = this.ReducedModel.System.f.M2 * v2;
%                 a_4 = 0;
%                 if ~isempty(fs.B)
%                     a_4 = (I-V*V')*fs.B.evaluate(t,mu)*ut;
%                 end
%                 a2 = sqrt((a_1+a_2-a_3+a_4)'*this.ReducedModel.FullModel.G...
%                     *(a_1+a_2-a_3+a_4));
%                 if abs((a-a2)/a2) > 1e-1
%                     error('alpha computation mismatch.');
%                 end
            end
        end
        
        function b = getBeta(this, x, t, mu)
            x = x(1:end-1);
            % Old log lip const kernel learning code
            %x = x .* (this.scale(:,2) - this.scale(:,1)) + this.scale(:,1);
            %             x = (x - this.scale(:,1)) ./ (this.scale(:,2) - this.scale(:,1));
            %             b = this.kexp.evaluate(x, t, mu)*olderr + sqrt(abs(a));
            
            % Validation code
            if this.UseTrueLogLipConst || this.UseJacobianLogLipConst
                rm = this.ReducedModel;
                f = rm.FullModel.System.f;
                tx = this.fullsol(:,find(rm.scaledTimes - t < eps,1));
                rx = rm.V*x;
                d = tx - rx;
                diff = sum(d.*d);
                if diff ~= 0
                    % Use local jacobian matrix log lip const
                    if this.UseJacobianLogLipConst
                        b = d'*(f.getStateJacobian(tx,t,mu)*d) / diff;
                    % Use local function log lip const
                    else
                        b = d'*(f.evaluate(tx,t,mu) - f.evaluate(rx,t,mu)) / diff;
                    end
                else
                    b = 0;
                end
            elseif this.UseFullJacobian
                rm = this.ReducedModel;
                f = rm.FullModel.System.f;
                b = general.Utils.logNorm(f.getStateJacobian(rm.V*x, t, mu));
            else
                DJ = this.JacMDEIM.evaluate(x,t,mu);
                b = general.Utils.logNorm(DJ);
            end
            
            %             DFJ = fm.Approx.getStateJacobian(rx,t,mu);
            %             jtlc = (d'*(FJ*d)) / sum(d.*d);
            %             djtlc = (d'*(DFJ*d)) / sum(d.*d);
            
            %             fprintf('True: %g, Jac: %g, diff: %g, rel: %g, DEIMJac: %g, diff: %g, rel: %g\n',...
            %                 tlc,jtlc,abs(tlc-jtlc),abs((tlc-jtlc)/tlc),...
            %                 djtlc,abs(tlc-djtlc),abs((tlc-djtlc)/tlc));
            
            %             J = this.deim.getStateJacobian(x,t,mu);
            %             jln = general.Utils.logNorm(J);
            
            %             fprintf('Log norms - F: %g, SF: %g, SR: %g, R: %g\n',...
            %                 fjln,sfjln,dsjln,jln); %, FJN: %g
            
            if isnan(b)
                fprintf(2,'NaN b value occured in error.DEIMEstimator.getBeta!\n');
                b = 0;
            end
            
            % Take care of the A(t,\mu) part, if existing
            %this.kexp(1,end+1) = b;
            if ~isempty(this.Aln)
                b = b + this.Aln.compose(t, mu);
            end
            %this.kexp(2,end) = b;
        end
        
        function eint = evalODEPart(this, x, t, mu, ut)
            % Compose the ode function
            eint = this.getBeta(x, t, mu)*x(end) + this.getAlpha(x, t, mu, ut);
        end
        
        function copy = clone(this)
            copy = clone@error.BaseEstimator(this, error.DEIMEstimator);
            copy.silent = true;
            
            % DEIM stuff
            copy.JacMatDEIMMaxOrder = this.JacMatDEIMMaxOrder;
            copy.JacMDEIM = [];
            if ~isempty(this.JacMDEIM)
                copy.JacMDEIM = this.JacMDEIM.clone;
            end
            
            % Sim Trans stuff
            copy.JacSimTransMaxSize = this.JacSimTransMaxSize;
            copy.jstSize = this.jstSize;
            copy.QFull = this.QFull;
            
            copy.UseTrueLogLipConst = this.UseTrueLogLipConst;
            copy.UseTrueDEIMErr = this.UseTrueDEIMErr;
            copy.UseFullJacobian = this.UseFullJacobian;
            
            copy.TrainDataSelector = this.TrainDataSelector.clone;
            
            copy.scale = this.scale;
            copy.kexp = this.kexp;
            copy.fullsol = this.fullsol;
            
            copy.M3 = this.M3;
            copy.M4 = this.M4;
            copy.M5 = this.M5;
            if ~isempty(this.M6)
                copy.M6 = this.M6.clone;
            end
            if ~isempty(this.M7)
                copy.M7 = this.M7.clone;
            end
            if ~isempty(this.M8)
                copy.M8 = this.M8.clone;
            end
            if ~isempty(this.M9)
                copy.M9 = this.M9.clone;
            end
            if ~isempty(this.M10)
                copy.M10 = this.M10.clone;
            end
            if ~isempty(this.M11)
                copy.M11 = this.M11.clone;
            end
            if ~isempty(this.M12)
                copy.M12 = this.M12.clone;
            end
            if ~isempty(this.Ah)
                copy.Ah = this.Ah.clone;
                copy.Aln = this.Aln;
            end
            if ~isempty(this.Bh)
                copy.Bh = this.Bh.clone;
            end
            copy.silent = false;
        end
        
    end
    
    methods(Access=protected)
        function time = postprocess(this, x, varargin)
            % Postprocessing for the error estimate.
            %
            % Parameters:
            % x: The whole just computed trajectory @type matrix<double>
            % varargin: Any further arguments that are passed, but ignored
            % here as not necessary. @type varargin
            %
            % t, mu, inputidx
            % do nothing
            time = tic;
            this.StateError = x(end,:);
            time = toc(time);
        end
    end
    
    methods(Access=private)
        function compute_jacmdeim(this, fm)
            %% Large offline data
            jtd = data.ApproxTrainData.computeFrom(fm, ...
                general.JacCompEvalWrapper(fm.System.f), this.TrainDataSelector, false);
            md = fm.Data;
            md.JacobianTrainData = jtd;
            
            %% Matrix DEIM
            jd = general.MatrixDEIM;
            jd.MaxOrder = this.JacMatDEIMMaxOrder;
            jd.NumRows = fm.System.f.fDim;
            if KerMor.App.Verbose > 0
                fprintf('Computing matrix DEIM (MaxOrder=%d) of system jacobian...\n',...
                    this.JacMatDEIMMaxOrder);
            end
            jd.computeDEIM(general.JacCompEvalWrapper(fm.System.f), ...
                md.JacobianTrainData.fxi);
            % Project, as arguments that will be passed are in reduced
            % dimension
            this.JacMDEIM = jd;
            % Initialize to certain order
            this.JacMatDEIMOrder = ceil(this.JacMatDEIMMaxOrder/2);
        end
        
        function compute_simtrans(this, fm)
            if ~isempty(this.JacSimTransMaxSize)
                jtd = fm.Data.JacobianTrainData;

                d = fm.System.f.xDim;
                n = size(jtd.fxi,2);
                v = data.FileMatrix(d,n,fm.Data.DataDirectory,512*1024^2);
                ln = zeros(1,n);
                times = ln;
                pi = tools.ProcessIndicator('Computing Jacobian similarity transform data for %d jacobians',n,false,n);
                hassparse = ~isempty(fm.System.f.JSparsityPattern);
                if hassparse
                    [i,j] = find(fm.System.f.JSparsityPattern);
                end
                for nr = 1:n
                    if hassparse
                        J = sparse(i,j,jtd.fxi(:,nr),d,d);
                    else
                        J = reshape(jtd.fxi(:,nr),d,d);
                    end
                    t = tic;
                    [ln(nr), v(:,nr)] = general.Utils.logNorm(J);
                    times(nr) = toc(t);
                    pi.step;
                end
                pi.stop;

                jstd.VFull = v;
                jstd.LogNorms = ln;
                jstd.CompTimes = times;
                fm.Data.JacSimTransData = jstd;
                
                p = general.POD;
                p.Mode = 'abs';
                p.Value = this.JacSimTransMaxSize;
                [Q, this.QSingVals] = p.computePOD(jstd.VFull);
                if isa(Q,'data.FileMatrix')
                    this.QFull = Q.toMemoryMatrix;
                else
                    this.QFull = Q;
                end
                if size(this.QFull,2) < this.JacSimTransMaxSize
                    fprintf(2,'Only %d nonzero singular values in difference to %d desired ones. Setting JacSimTransMaxSize=%d\n',...
                        size(this.QFull,2),this.JacSimTransMaxSize,size(this.QFull,2));
                    this.silent = true;
                    this.JacSimTransMaxSize = size(this.QFull,2);
                    this.silent = false;
                end
                red = 1-size(this.QFull,2)/size(this.QFull,1);
                if KerMor.App.Verbose > 0
                    fprintf(['Computed partial similarity transform with target size %d over %d eigenvectors. '...
                        'Resulting reduction %d/%d (%g%%)\n'],...
                        p.Value,size(jstd.VFull,2),...
                        size(this.QFull,2),size(this.QFull,1),100*red);
                end
                if red < .2
                    this.JacSimTransSize = [];
                    fprintf(2,'Achieved reduction under 20%%, disabling similarity transformation.\n');
                else
                    this.JacSimTransSize = ceil(this.JacSimTransMaxSize/2);
                end
            end
        end
        
        function handleOrderUpdate(this, ~, ~)
            this.updateErrMatrices;
        end
        
        function updateErrMatrices(this)
            % This methods updates all the offline-computable matrices
            % involved for error estimation that depend on the DEIM
            % approximation order.
            if isempty(this.deim)
                if isempty(this.ReducedModel)
                    error('Cannot update error matrices as no reduced model instance is specified for this estimator.');
                else
                    error('Unintended error. No deim instance stored but a reduced model is set.');
                end
            end
            
            G = this.ReducedModel.FullModel.G;
            
            if KerMor.App.Verbose > 3
                fprintf('error.DEIMEstimator: Updating error matrices of DEIMEstimator (#%s) to [%d %d]\n',...
                    this.ID,this.deim.Order);
            end
            
            if this.deim.Order(2) > 0
                M1 = this.deim.M1;
                M2 = this.deim.M2;
                this.M3 = M1'*(G*M1);
                this.M4 = 2*M1'*(G*M2);
                this.M5 = M2'*(G*M2);
                if ~isempty(this.Ah)
                    this.M7 = 2*M1'*(G*this.Ah);
                    this.M8 = 2*M2'*(G*this.Ah);
                end
                if ~isempty(this.Bh)
                    this.M10 = 2*M1'*(G*this.Bh);
                    this.M11 = 2*M2'*(G*this.Bh);
                end
            else
                this.M3 = [];
                this.M4 = [];
                this.M5 = [];
                this.M7 = [];
                this.M8 = [];
                this.M10 = [];
                this.M11 = [];
            end
        end
        
        function computeLogNormApprox(this)
            %             [res, mScale, MScale] = testing.LogNorm(this.ReducedModel.FullModel);
            %             a = approx.algorithms.VectorialKernelOMP;
            %             a.NumGammas = 60;
            %             a.MinGFactor = .5;
            %             a.MaxGFactor = 7;
            %             a.UseOGA = true;
            %             a.UsefScaling = true;
            %             k = kernels.KernelExpansion;
            %             k.Kernel = kernels.GaussKernel;
            %             a.computeApproximation(k, res);
            %             this.kexp = k;
            %             this.scale = [mScale MScale];
            
            load(fullfile(KerMor.App.HomeDirectory,'data/lognorm/res_loclognorms_scaled.mat'));
            this.kexp = kexp;%#ok
            this.scale = [ms Ms];
        end
    end
    
    %% Getter & Setter
    methods
        function value = get.JacSimTransSize(this)
            value = this.jstSize;
        end
        
        function set.JacSimTransSize(this, value)
            if isempty(value) || (isposintscalar(value) ...
                    && value <= this.JacSimTransMaxSize)
                if ~isequal(this.jstSize,value)
                    this.jstSize = value;
                    this.JacMDEIM.setSimilarityTransform(...
                        this.QFull(:,1:this.jstSize));
                end
            else
                error('Value must be empty or a positive int scalar smaller than JacSimTransMaxSize');
            end
        end
        
        function value = get.JacMatDEIMOrder(this)
            if ~isempty(this.JacMDEIM)
                value = this.JacMDEIM.Order(1);
            else
                value = -1;
            end
        end
        
        function set.JacMatDEIMOrder(this, value)
            if ~isposintscalar(value)
                error('JacMatDEIMOrder must be a positive integer.');
            end
            if value <= this.JacMatDEIMMaxOrder
                if this.JacMDEIM.Order(1) ~= value
                    this.JacMDEIM.Order = value;
                end
            else
                error('Cannot set an order higher (%d) than the current JacMatDEIMMaxOrder of %d',...
                    value,this.JacMatDEIMMaxOrder);
            end
        end
        
        function set.TrainDataSelector(this, value)
            this.checkType(value, 'data.selection.ASelector');
            this.TrainDataSelector = value;
        end
        
        function set.JacSimTransMaxSize(this, value)
            if ~this.silent && ~isempty(value) && ~isempty(this.JacSimTransMaxSize) && ...
                    this.JacSimTransMaxSize ~= value %#ok<MCSUP>
                fprintf(2,'New maximum size of similarity transform. Re-run of offlineComputations necessary.\n');
            end
            this.JacSimTransMaxSize = value;
        end
    end
    
     methods(Static, Access=protected)
        function obj = loadobj(obj)
            if ~isempty(obj.deim)
                obj.uolst = addlistener(obj.deim,'OrderUpdated',@obj.handleOrderUpdate);
            end
        end
     end
    
    methods(Static)
        
        % Abstract static method that forces subclasses to specify whether
        % an estimator can be used for a given model or not.
        function errmsg = validModelForEstimator(model)
            errmsg = '';
            if ~isa(model.Approx,'approx.DEIM')
                errmsg = 'The reduced models core function must be a DEIM approximation.';
            elseif ~isempty(model.System.B) && ...
                    ~any(strcmp(class(model.System.B),{'dscomponents.LinearInputConv','dscomponents.AffLinInputConv'}))
                errmsg = 'The systems input must be either Linear- or AffLinInputConv';
            elseif ~isempty(model.System.A) && ...
                    ~any(strcmp(class(model.System.A),...
                    {'dscomponents.LinearCoreFun','dscomponents.AffLinCoreFun'}))
                errmsg = ['If the system component A is set, is has to be either a'...
                    'LinearCoreFun or AffLinCoreFun.'];
            end
        end
    end
end