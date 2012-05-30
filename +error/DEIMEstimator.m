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
        
        % Stores the projected DEIM approximation
        deim;
        
        % Stores zi min- and max scaling
        scale;
        
        % The local logarithmic norm estimation
        kexp;
    end
    
    methods
        function this = DEIMEstimator(rmodel)
            this = this@error.BaseEstimator;
            this.ExtraODEDims = 1;
            if nargin == 1    
                this.setReducedModel(rmodel);
            end
        end
        
        function setReducedModel(this, rmodel)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            
            % Call superclass
            setReducedModel@error.BaseEstimator(this, rmodel);
            fm = rmodel.FullModel;
            this.deim = rmodel.System.f;
            A = [];
            if ~isempty(fm.System.A)
                A = (fm.System.A - rmodel.V*(rmodel.W'*fm.System.A))*rmodel.V;
                this.M6 = A'*A;
                this.Ah = A;
            else
                this.M6 = []; this.M7 = []; this.M8 = []; this.M9 = [];
            end
            B = [];
            if ~isempty(fm.System.B)
                B = fm.System.B - rmodel.V*(rmodel.W'*fm.System.B);
                this.M12 = B'*B;
                this.Bh = B;
            else
                this.M9 = []; this.M10 = []; this.M11 = []; this.M12 = [];
            end
            if ~isempty(A) && ~isempty(B)
                this.M9 = 2*A'*B;
            else
                this.M9 = [];
            end
            
            % Update the other matrices and add a listener to the deim
            % instance to react to DEIM order changes
            this.updateErrMatrices;    
            addlistener(this.deim,'Order','PostSet',@this.orderUpdated);
            
            % Compute approximation for local logarithmic norm
%             this.computeLogNormApprox;
        end
        
        function eint = evalODEPart(this, x, t, mu, ut)
            olderr = x(end);
            x = x(1:end-1);
            v1 = this.deim.f.evaluateComponentSet(1, x, t, mu);
            v2 = this.deim.f.evaluateComponentSet(2, x, t, mu);
            
            % NOTE: The factors 2 before some summands are added at offline
            % stage!!
            a = v1'*this.M3*v1 - v1'*this.M4*v2 + v2'*this.M5*v2;
            if ~isempty(this.M6) % An time/param affine A is set
                a = a + x'*this.M6.evaluate(t,mu)*x + x'*this.M7.evaluate(t,mu)*v1 ...
                    - x'*this.M7.evaluate(t,mu)*v2;
                if ~isempty(ut)
                    a = a + x'*this.M9.evaluate(t,mu)*ut;
                end
            end
            if ~isempty(ut) % An input function u is set
                a = a + v1'*this.M10.evaluate(t,mu)*ut - v2'*this.M11.evaluate(t,mu)*ut ...
                    + ut'*this.M12.evaluate(t,mu)*ut;
            end
            %x = x .* (this.scale(:,2) - this.scale(:,1)) + this.scale(:,1);
%             x = (x - this.scale(:,1)) ./ (this.scale(:,2) - this.scale(:,1));
%             eint = this.kexp.evaluate(x, t, mu)*olderr + sqrt(abs(a));
            J = this.deim.getStateJacobian(x,t,mu);
            b = general.Utils.logNorm(J);
            eint = b*olderr + sqrt(abs(a));
        end
        
        function copy = clone(this)
            copy = clone@error.BaseEstimator(this, error.DEIMEstimator(this.ReducedModel));
            copy.Uerr1 = this.Uerr1;
            copy.deim = this.deim.clone;
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
        
        function orderUpdated(this, ~, ~)%sender, data
            if KerMor.App.Verbose > 0
                fprintf('Updating DEIM error estimator matrices..\n');
            end
            if ~this.Enabled && this.deim.Order(2) > 0
                this.Enabled = true;
                fprintf('DEIM error estimation order set. Enabling DEIMEstimator.\n');
            end
            this.updateErrMatrices;
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
    
    methods(Static)
        % Abstract static method that forces subclasses to specify whether
        % an estimator can be used for a given model or not.
        function errmsg = validModelForEstimator(rmodel)
            errmsg = '';
            if ~isa(rmodel.System.f,'approx.DEIM')
                errmsg = 'The reduced models core function must be a DEIM approximation.';
            elseif ~isempty(rmodel.FullModel.System.B) && ...
                   ~any(strcmp(class(rmodel.FullModel.System.B),{'dscomponents.LinearInputConv','dscomponents.AffLinInputConv'}))
                errmsg = 'The systems input must be either Linear- or AffLinInputConv';
            end
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj)
            obj = loadobj@error.BaseEstimator(obj);
            addlistener(obj.deim,'Order','PostSet',@obj.orderUpdated);
        end
    end
end