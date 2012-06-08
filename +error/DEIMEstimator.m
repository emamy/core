classdef DEIMEstimator < error.BaseEstimator
% DEIMEstimator: A-posteriori error estimation for DEIM reduced models.
%
% @author Daniel Wirtz @date 2012-05-10
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
        % @type integer @default []
        %
        % See also: JacSimTransSize
        JacSimTransMaxSize = [];
        
        % The selector that chooses the full model's trajectory points that
        % are used for the MatrixDEIM approximation of the jacobian and the
        % partial similarity transformation.
        %
        % @type approx.selection.ASelector @default
        % approx.selection.DefaultSelector
        TrainDataSelector;
    end
    
    properties(SetAccess=private)
        % The MatrixDEIM instance used to approximate the state space
        % jacobian matrices
        %
        % The default MaxOrder setting is 50.
        %
        % @type general.MatrixDEIM
        JacMDEIM;
        
        % The POD instance to compute the partial similarity transform
        % matrix.
        %
        % The POD settings default to Mode 'eps' with Value = 1e-7, causing
        % the JacSimTransMaxSize to be set to the resulting POD size.
        % If the Mode is set to 'abs', the current value of
        % JacSimTransMaxSize will be used in offlineGenerations.
        %
        % @type general.POD
        SimTransPOD;
                
        % The complete similarity transformation matrix of size 
        % `d \times JacSimTransMaxSize`
        %
        % (Debug/testing use)
        %
        % @type matrix<double> @default []
        STFull = [];
        
        % The singular values computed by the SimTransPOD in the offline
        % stage.
        %
        % (Debug/testing use)
        %
        % @type rowvec<double> @default []
        STSingVals = [];
        
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
        % The DEIM error estimation order to use for simulations.
        %
        % Setting to a positive value results in activation of the error
        % estimator.
        %
        % @type integer @default 5
        ErrorOrder;
        
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
    end
    
    properties(Access=private)
        eOrder = 5;
        jstOrder = 0;
    end
    
    properties(Access=private, Transient)
        fullsol;
        
        % Stores the reduced system's deim instance
        deim;
    end
    
    methods
        function this = DEIMEstimator
            this = this@error.BaseEstimator;
            this.ExtraODEDims = 1;
            this.JacMDEIM = general.MatrixDEIM;
            this.JacMDEIM.MaxOrder = 50;
            p = general.POD;
            p.Mode = 'eps';
            p.Value = 1e-7;
            this.SimTransPOD = p;
            this.TrainDataSelector = approx.selection.DefaultSelector;
        end
        
        function offlineComputations(this, fm)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            
            % Call superclass
            offlineComputations@error.BaseEstimator(this, fm);
            
            fA = fm.System.A;
            if fm.GScaled ~= 1
                error('Not yet implemented for non-1 G matrices.');
            end
            md = fm.Data;
            A = [];
            if ~isempty(fA)
                A = (fA - md.V*(md.W'*fA))*md.V;
                this.M6 = A'*A;
                this.Ah = A;
                
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
            else
                this.M6 = []; this.M7 = []; this.M8 = []; this.M9 = [];
            end
            B = [];
            if ~isempty(fm.System.B)
                B = fm.System.B - md.V*(md.W'*fm.System.B);
                this.M12 = B'*B;
                this.Bh = B;
            else
                this.M9 = []; this.M10 = []; this.M11 = []; this.M12 = [];
            end
            if ~isempty(A) && ~isempty(B)
                if isa(B,'dscomponents.LinearInputConv')
                    this.M9 = 2*A'*B.B;
                else
                    % Both inherit from AffParamMatrix, so they can be
                    % multiplied :-)
                    this.M9 = 2*A'*B;
                end
            else
                this.M9 = [];
            end
            
            % LogNorm-related computations
            jd = this.JacMDEIM;
            jd.NumRows = fm.System.f.XDim;
            if KerMor.App.Verbose > 0
                fprintf('Computing matrix DEIM (MaxOrder=%d) of system jacobian...\n',...
                    this.JacMDEIM.MaxOrder);
            end
            jd.computeDEIM(general.JacCompEvalWrapper(fm.System.f), ...
                md.JacobianTrainData.fxi);
            % Project, as arguments that will be passed are in reduced
            % dimension
            this.JacMDEIM = jd.project(md.V, md.W);
            
            if ~isempty(md.JacSimTransData.VFull)
                p = this.SimTransPOD;
                if strcmp(p.Mode,'abs')
                    if isempty(this.JacSimTransMaxSize)
                        error('When SimTransPOD has mode ''abs'' the property JacSimTransMaxSize must be set.');
                    end
                    p.Value = this.JacSimTransMaxSize;
                end
                if KerMor.App.Verbose > 0
                    fprintf(['Computing partial similarity transform with POD Mode=''%s'' '...
                        'and Value=%g on %d eigenvectors...\n'],...
                        p.Mode,p.Value,size(md.JacSimTransData.VFull,2));
                end
                [this.STFull, this.STSingVals] = ...
                    p.computePOD(md.JacSimTransData.VFull);
                this.JacSimTransMaxSize = size(this.STFull,2);
                this.JacSimTransSize = ceil(this.JacSimTransMaxSize/2);
            else
                fprintf(2,'ModelData.JacSimTransData.VFull is empty. Not computing similarity transform.\n');
            end
            
            % Compute approximation for local logarithmic norm
%             this.computeLogNormApprox;
        end
        
        function ct = prepareConstants(this, mu, inputidx)
            % Experimental part
            [~, this.fullsol] = this.ReducedModel.FullModel.computeTrajectory(mu,inputidx);
            
            t = tic;
            if isempty(this.deim)
                this.deim = this.ReducedModel.System.f;
            end
            neworder = [this.deim.Order(1) this.ErrorOrder];
            % Only update DEIM error matrices if a new order has been set.
            if ~isequal(this.deim.Order, neworder)
                this.deim.Order = neworder;
                this.updateErrMatrices;
            end
            ct = toc(t);
        end
        
        function a = getAlpha(this, x, t, mu, ut)
            x = x(1:end-1);
            v1 = this.deim.f.evaluateComponentSet(1, x, t, mu);
            v2 = this.deim.f.evaluateComponentSet(2, x, t, mu);
            
            % NOTE: The factors 2 before some summands are added at offline
            % stage!!
            a = v1'*this.M3*v1 - v1'*this.M4*v2 + v2'*this.M5*v2;
%             fprintf('a: %g, a_1:%g, a_2:%g, a_3:%g, M3: %dx%d, M4: %dx%d, M5:%dx%d\n',a,...
%                 v1'*this.M3*v1, - v1'*this.M4*v2, v2'*this.M5*v2,...
%                 size(this.M3,1),size(this.M3,2),size(this.M4,1),size(this.M4,2),...
%                 size(this.M5,1),size(this.M5,2));
            
            if ~isempty(this.M6) % An time/param affine A is set
                a = a + x'*this.M6.evaluate(x,t,mu) + v1'*this.M7.evaluate(x,t,mu) ...
                    - v2'*this.M8.evaluate(x,t,mu);
%                 fprintf('a+A-pt: %g, M6:%g, M7:%g, M8:%g\n',a,...
%                     x'*this.M6.evaluate(x,t,mu), v1'*this.M7.evaluate(x,t,mu), ...
%                     - v2'*this.M8.evaluate(x,t,mu));
                if ~isempty(ut)
                    a = a + x'*this.M9.evaluate(ut,t,mu);
                end
            end
            if ~isempty(ut) % An input function u is set
                a = a + v1'*this.M10.evaluate(t,mu)*ut - v2'*this.M11.evaluate(t,mu)*ut ...
                    + ut'*this.M12.evaluate(t,mu)*ut;
            end
            a = sqrt(abs(a));
            
            %% Validation: direct computation (expensive)
%             V = this.ReducedModel.V;
%             fs = this.ReducedModel.FullModel.System;
%             I = speye(size(V,1));
%             a1 = (I-V*V')*fs.A.evaluate(V*x,t,mu);
%             a2 = this.deim.M1 * v1;
%             a3 = this.deim.M2 * v2;
%             a2 = norm(a1+a2-a3);
        end
        
        function b = getBeta(this, x, t, mu)
            x = x(1:end-1);
            %x = x .* (this.scale(:,2) - this.scale(:,1)) + this.scale(:,1);
%             x = (x - this.scale(:,1)) ./ (this.scale(:,2) - this.scale(:,1));
%             eint = this.kexp.evaluate(x, t, mu)*olderr + sqrt(abs(a));
            
            % Validation code
            rm = this.ReducedModel;
            fm = rm.FullModel;
            ff = fm.System.f;
%             tx = this.fullsol(:,find(rm.scaledTimes - t < eps,1));
            rx = rm.V*x;
%             d = tx - rx;
%             tlc = d'*(ff.evaluate(tx,t,mu) - ff.evaluate(rx,t,mu)) ...
%                 / sum(d.*d);
            
            FJ = ff.getStateJacobian(rx,t,mu);
%             DFJ = fm.Approx.getStateJacobian(rx,t,mu);
%             jtlc = (d'*(FJ*d)) / sum(d.*d);
%             djtlc = (d'*(DFJ*d)) / sum(d.*d);
            
%             fprintf('True: %g, Jac: %g, diff: %g, rel: %g, DEIMJac: %g, diff: %g, rel: %g\n',...
%                 tlc,jtlc,abs(tlc-jtlc),abs((tlc-jtlc)/tlc),...
%                 djtlc,abs(tlc-djtlc),abs((tlc-djtlc)/tlc));
            
            fjln = general.Utils.logNorm(FJ);
            ST = this.JacMDEIM.ST;
            sfjln = general.Utils.logNorm(ST'*(FJ*ST));
            
            DJ = this.JacMDEIM.evaluate(x,t,mu);
            dsjln = general.Utils.logNorm(DJ);
            
            J = this.deim.getStateJacobian(x,t,mu);
            jln = general.Utils.logNorm(J);
            
            fprintf('Log norms - F: %g, SF: %g, SR: %g, R: %g\n',...
                fjln,sfjln,dsjln,jln); %, FJN: %g

%             b = tlc;
            b = fjln;
            
            if isnan(b)
                b = 0;
            end
            
            % Take care of the A(t,\mu) part, if existing
            if ~isempty(this.Aln)
                b = b + this.Aln.compose(t, mu)*this.ReducedModel.dt;
            end
        end
        
        function eint = evalODEPart(this, x, t, mu, ut)
            % Compose the ode function
            eint = this.getBeta(x, t, mu)*x(end) + this.getAlpha(x, t, mu, ut);
        end
        
        function copy = clone(this)
            copy = clone@error.BaseEstimator(this, error.DEIMEstimator(this.ReducedModel));
            
            % Dont copy the DEIM approx function, the're working on the same!
            copy.deim = this.deim;
            
            copy.JacMDEIM = this.JacMDEIM.clone;
            copy.JacSimTransMaxSize = this.JacSimTransMaxSize;
            copy.eOrder = this.eOrder;
            copy.jstOrder = this.jstOrder;
            copy.STFull = this.STFull;
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
        end
        
        function value = get.ErrorOrder(this)
            value = this.eOrder;
        end
        
        function set.ErrorOrder(this, value)
            if value == 0
                this.Enabled = false;
            elseif ~isposintscalar(value)
                error('The error order has to be an integer scalar.');
            else
                this.Enabled = true;
            end
            this.eOrder = value;
        end
        
        function value = get.JacSimTransSize(this)
            value = this.jstOrder;
        end
        
        function set.JacSimTransSize(this, value)
            if value < this.JacSimTransMaxSize
                if this.jstOrder ~= value
                    this.jstOrder = value;
                    this.JacMDEIM.setSimilarityTransform(...
                        this.STFull(:,1:this.jstOrder));
                end
            end
        end
        
        function set.TrainDataSelector(this, value)
            this.checkType(value, 'approx.selection.ASelector');
            this.TrainDataSelector = value;
        end
        
        function set.JacSimTransMaxSize(this, value)
            if ~isempty(this.JacSimTransMaxSize) && ...
                    this.JacSimTransMaxSize ~= value
                fprintf('New maximum size of similarity transform. Re-run of offlineComputations necessary.\n');
            end
            this.JacSimTransMaxSize = value;
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
        function updateErrMatrices(this)
            % This methods updates all the offline-computable matrices
            % involved for error estimation that depend on the DEIM
            % approximation order.
            if this.deim.Order(2) > 0
                M1 = this.deim.M1;
                M2 = this.deim.M2;
                this.M3 = M1'*M1;
                this.M4 = 2*M1'*M2;
                this.M5 = M2'*M2;
                if ~isempty(this.Ah)
                    this.M7 = 2*M1'*this.Ah;
                    this.M8 = 2*M2'*this.Ah;
                end
                if ~isempty(this.Bh)
                    this.M10 = 2*M1'*this.Bh;
                    this.M11 = 2*M2'*this.Bh;
                end
            else
                if this.Enabled
                    this.Enabled = false;
                    fprintf('DEIM error estimation order not set. Disabling DEIMEstimator.\n');
                end
            end
        end
        
%         function orderUpdated(this, ~, ~)%sender, data
%             if KerMor.App.Verbose > 0
%                 fprintf('Updating DEIM error estimator matrices.. (ID %s)\n',this.ID);
%             end
%             if ~this.Enabled && this.deim.Order(2) > 0
%                 this.Enabled = true;
%                 fprintf('DEIM error estimation order set. Enabling DEIMEstimator.\n');
%             end
%             this.updateErrMatrices;
%         end
        
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