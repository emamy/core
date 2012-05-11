classdef DEIMEstimator < error.BaseEstimator
% DEIMEstimator: 
%
%
%
% @author Daniel Wirtz @date 2012-05-10
%
% @new{0,6,dw,2012-05-10} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        M1 = [];
        M2 = [];
        M3 = [];
        M4 = [];
        B = [];
        
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
            this.deim = rmodel.System.f;
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
            if ~isempty(fm.System.B)
                this.B = fm.System.B - rmodel.V*(rmodel.W'*fm.System.B);
                this.M4 = this.B'*this.B;
            end
            
            addlistener(rmodel.System.f,'Order','PostSet',@this.orderUpdated);
            this.updateErrMatrices;
            
            this.computeLogNormApprox;
        end
        
        function eint = evalODEPart(this, x, t, mu, ut)
            olderr = x(end);
            x = x(1:end-1);
            v1 = this.deim.f.evaluateComponentSet(1, x, t, mu);
            v2 = this.deim.f.evaluateComponentSet(2, x, t, mu);
            a = v1'*v1 + v2'*v2 - 2*v1'*this.M1*v2;
            if ~isempty(ut) % An input function u is set
                a = a + v1'*this.M2.evaluate(t,mu)*ut - v2'*this.M3.evaluate(t,mu)*ut ...
                    + ut'*this.M4.evaluate(t,mu)*ut;
            end
            x = x .* (this.scale(:,2) - this.scale(:,1)) + this.scale(:,1);
            eint = this.kexp.evaluate(x, t, mu)*olderr + sqrt(abs(a));
        end
        
        function copy = clone(this)
            copy = clone@error.BaseEstimator(this, error.DEIMEstimator(this.ReducedModel));
            copy.Uerr1 = this.Uerr1;
            copy.deim = this.deim.clone;
        end
    end
    
    methods(Access=protected)
        function time = postprocess(this, x, varargin)
            % t, mu, inputidx
            % do nothing
            time = tic;
            this.StateError = x(end,:);
            time = toc(time);
        end
    end
    
    methods(Access=private)
        function updateErrMatrices(this)
            d = this.deim;
            this.M1 = d.Uerr1'*d.Uerr2;
            if ~isempty(this.B)
                this.M2 = d.Uerr1'*this.B;
                this.M3 = d.Uerr2'*this.B;
            end
        end
        
        function orderUpdated(this, ~, ~)%sender, data
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
            obj = loadobj@BaseEstimator(obj);
            addlistener(obj.ReducedModel.System.f,'Order','PostSet',@obj.orderUpdated);
        end
    end
end