classdef TPWLApprox < approx.BaseApprox
% Trajectory-piecewise function approximation
%
% Implements the trajectory-piecewise approximation algorithm proposed in [RW01/Re03],
% Rewienski, M.: A trajectory piecewise-linear approach to model order reduction of nonlinear
% dynamical systems. Ph.D. thesis, Citeseer (2003)
%
% The implementation is only completely like the proposed one if the data.selection.EpsSelector is
% used in conjunction with this class.
%
% @author Daniel Wirtz @date 2011-04-01
%
% @change{0,5,dw,2011-07-07} Changed this class to adopt to the new BaseApprox interface.
%
% @change{0,4,dw,2011-05-10}
% - The `A_i` matrices are now sparse if applicable (numel ~= 0 / numel < .5)
% - Removed the WeightFun property and reintroduced the Beta property, correctly forwarding the
% setting to the internally used gauss kernel.
%
% @new{0,4,dw,2011-05-06} Finished the TPWL implementation. Multiargument-evaluations now work
% and the WeightFun was set to the TPWL source papers default `e^{-\beta \frac{d_i}{m}}`.
%
% @change{0,3,sa,2011-04-21} Implemented Setters for the properties
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetObservable)        
        % The minimum value a weight needs to have after normalization in
        % order to affect the evaluation. Set this value to a certain
        % threshold in order to restrict globally supported weight
        % functions (i.e. the Gaussian) to a local support.
        %
        % @propclass{experimental} Might not be used in the end, so far only for speedup (locality)
        MinWeightValue = 1e-10;
    end
    
    properties(Dependent, SetObservable)
        % @propclass{critical}
        Beta;
    end
    
    properties(SetAccess=private)
        % The centers
        xi;
        % The gradient matrix
        Ai;
        bi;
        GaussWeight;
        fBeta = 25;
    end
    
    methods        
        function set.MinWeightValue(this, value)
            if ~isposrealscalar(value)
                error('MinWeightValue must be a positive real scalar');
            end
            this.MinWeightValue = value;
        end
        
        function set.Beta(this, value)
            if ~isposrealscalar(value)
                error('Beta must be a positive real scalar.');
            end
            this.fBeta = value;
            this.GaussWeight.Gamma = 1/sqrt(value);
        end
        
        function value = get.Beta(this)
            value = this.fBeta;
        end
                
        function this = TPWLApprox(sys)
            this = this@approx.BaseApprox(sys);
            
            % Beta from Paper is e^{-25 d/m}, for gauss: e^{d/(m*g^2)}, so g=.2
            % Initialize the gaussian weight with what ever Beta is set to
            this.GaussWeight = kernels.GaussKernel(1/sqrt(this.Beta));
            this.registerProps('Beta','MinWeightValue');
            
            % Set training data selection to epsselector (default strategy for original TPWL)
            this.TrainDataSelector = data.selection.EpsSelector;
            this.CustomProjection = true;
        end
        
        function copy = clone(this)
            copy = approx.TPWLApprox(this.System);
            copy = clone@approx.BaseApprox(this, copy);
            copy.Ai = this.Ai;
            copy.bi = this.bi;
            copy.xi = this.xi;
            copy.Beta = this.Beta;
            copy.MinWeightValue = this.MinWeightValue;
        end
        
        function projected = project(this, V, W)
            % Implements AProjectable.project()
            projected = project@approx.BaseApprox(this, V, W, this.clone);
            for i=1:length(this.Ai)
                projected.Ai{i} = W'*this.Ai{i}*V;
            end
            projected.bi = W'*this.bi;
            projected.xi = W'*this.xi;
        end
        
        function y = evaluate(this, x, t, mu)%#ok
            % Directly overrides the evaluate method of the ACoreFun as custom projection and
            % multi-argument evaluations are natively possible.
            % 
            % Implements ACoreFun abstract template method
            as = length(this.Ai);
            y = zeros(size(x));
            % Single argument evaluation
            if (size(x,2) == 1)
                % Compute all distances
                di = sqrt(sum((this.xi - repmat(x,1,as)).^2));
                di(di == 0) = eps;
                w = this.GaussWeight.evaluateScalar(sqrt(di/min(di)));
                % Normalize local weights
                w = w / sum(w);
                idx = 1:as;
                idx = idx(w > this.MinWeightValue);
                for i=idx
                    y = y + w(i)*(this.Ai{i}*x + this.bi(:,i));
                end
            % Multi argument evaluation
            else
                n = size(x,2);
                xisel = reshape(meshgrid(1:n,1:as),[],1)';
                % Compute all distances
                di = sqrt(sum((repmat(this.xi,1,n)-x(:,xisel)).^2));
                di = reshape(di',as,[])';
                di(di == 0) = eps;
                dimin = min(di,[],2);
                w = this.GaussWeight.evaluateScalar(sqrt(di./repmat(dimin,1,as)));
                w = w ./ repmat(sum(w,2),1,as);
                y = zeros(size(x));
                for i=1:as
                    % Find relevant weights
                    rel = w(:,i) > this.MinWeightValue;
                    % Normalize local weights
                    %wl = w(rel,i) ./ s(rel);
                    % Only update entries for relevant weights
                    hlp = this.Ai{i}*x(:,rel) + this.bi(:,i*ones(1,sum(rel)));
                    y(:,rel) = y(:,rel) + hlp*diag(w(rel,i));
                end
            end
        end
        
        function approximateSystemFunction(this, model)
            atd = model.Data.ApproxTrainData;
            this.xi = atd.xi;
            ti = atd.ti;
            mui = atd.mui;
            fxi = atd.fxi;
            
            as = size(this.xi,2);
            this.Ai = cell(0,as);
            for i = 1:as
                if isempty(mui)
                    mu = [];
                else
                    mu = mui(:,i);
                end
                model.System.setConfig(mu,model.DefaultInput);
                this.Ai{i} = model.System.f.getStateJacobian(this.xi(:,i),ti(i));
                this.bi(:,i) = fxi(:,i) - this.Ai{i}*this.xi(:,i);
            end
        end
        
        function y = evaluateCoreFun(this)%#ok
            % Noting to do here, evaluate is implemented directly. This method will never be called.
        end
    end
    
    methods(Static)
        function res = test_TWPLApprox
            m = models.synth.KernelTest(100,false);
            m.ErrorEstimator = [];
            m.Approx = approx.TPWLApprox(m.System);
            m.Approx.Beta = 25;
            m.dt = 0.01;
            m.Sampler = sampling.ManualSampler(m.getRandomParam);
            
            m.offlineGenerations;
            mu = m.getRandomParam;
            [~,y] = m.simulate(mu);
            r = m.buildReducedModel;
            [~,yr] = r.simulate(mu);
            fprintf('Max absolute error: %g\n',max(Norm.L2(y-yr)));
            res = true;
        end
    end
    
end

