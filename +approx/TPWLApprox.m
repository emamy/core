classdef TPWLApprox < approx.BaseApprox
% Trajectory-piecewise function approximation
%
% Implements the trajectory-piecewise approximation algorithm proposed in [RW01/Re03],
% Rewienski, M.: A trajectory piecewise-linear approach to model order reduction of nonlinear
% dynamical systems. Ph.D. thesis, Citeseer (2003)
%
% The implementation is only completely like the proposed one if the approx.selection.EpsSelector is
% used in conjunction with this class.
%
% @author Daniel Wirtz @date 2011-04-01
%
% @change{0,4,dw,2011-05-10}
% - The `A_i` matrices are now sparse if applicable (numel ~= 0 / numel < .5)
% - Removed the WeightFun property and reintroduced the Beta property, correctly forwarding the
% setting to the internally used gauss kernel.
%
% @new{0,4,dw,2011-05-06} Finished the TPWL implementation. Multiargument-evaluations now work
% and the WeightFun was set to the TPWL source papers default `e^{-\beta \frac{d_i}{m}}`.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
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
        Beta = 25;
    end
    
    properties(SetAccess=private)
        % The centers
        xi;
        % The gradient matrix
        Ai;
        bi;
        GaussWeight;
    end
    
    methods
        
        function this = TPWLApprox
            this = this@approx.BaseApprox;
            
            % Beta from Paper is e^{-25 d/m}, for gauss: e^{d/(m*g^2)}, so g=.2
            this.GaussWeight = kernels.GaussKernel(1/sqrt(this.Beta));
            this.registerProps('Beta','MinWeightValue');
            
            % Set training data selection to epsselector (default strategy for original TPWL)
            this.TrainDataSelector = approx.selection.EpsSelector;
            this.MultiArgumentEvaluations = true;
        end
        
        function copy = clone(this)
            % Implements ICloneable.clone()
            copy = approx.TPWLApprox;
            copy = clone@approx.BaseApprox(this, copy);
            copy.Ai = this.Ai;
            copy.bi = this.bi;
            copy.xi = this.xi;
            copy.Beta = this.Beta;
%             if isa(this.WeightFun,'ICloneable')
%                 copy.WeightFun = this.WeightFun.clone;
%             else
%                 warning('KerMor:cloning','Couldn''t properly clone the WeightFun as it''s not a ICloneable.');
%                 copy.WeightFun = this.WeightFun;
%             end
        end
        
        function projected = project(this, V, W)
            % Implements IProjectable.project()
            projected = this.clone;
            for i=1:length(this.Ai)
                projected.Ai{i} = W'*this.Ai{i}*V;
            end
            projected.bi = W'*this.bi;
            projected.xi = W'*this.xi;
        end
        
        function y = evaluateCoreFun(this, x, t, mu)
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
        
        function approximateCoreFun(this, model)
            % Implements BaseApprox abstract template method
            %
            % @todo create sparse matrices if suitable!
            
            % Load snapshots
            atd = model.Data.ApproxTrainData;
            
            % Compile necessary data
            this.xi = atd(4:end,:);
            ti = atd(3,:);
            muidx = atd(1,:);
            if all(muidx == 0)
                mui = [];
            else
                mui = model.Data.ParamSamples(:,muidx);
            end

            as = size(this.xi,2);
            N = size(this.xi,1);
            h = 1e-8;
            
            this.Ai = cell(0,as);
            dh = diag(ones(1,N)*h);
            for i = 1:as
                sel = ones(1,N)*i;
                xipt = this.xi(:,sel);
                if isempty(mui); mu = []; else mu = mui(:,sel); end
                tmp = (model.System.f.evaluate(xipt+dh,ti(sel),mu)-model.Data.ApproxfValues(:,sel))/h;
                % Make sparse if applicable
                if sum(sum(tmp ~= 0)) / numel(tmp) < .5
                    tmp = sparse(tmp);
                end
                this.Ai{i} = tmp;
                this.bi(:,i) = model.Data.ApproxfValues(:,i) - tmp*this.xi(:,i);
            end
        end
    end
    
    methods(Static)
        function res = test_TWPLApprox
            dims = 100;
            epsrad = .05;
            
            m = models.synth.KernelTest(dims,false);
            m.Approx = approx.TPWLApprox;
            %m.Approx.EpsRad = epsrad;
            d = sqrt(dims)*epsrad;
            k = kernels.GaussKernel;
            k.setGammaForDistance(50*d);
            m.Approx.WeightFun = k;
            m.dt = 0.01;
            
            m.offlineGenerations;
            %m.System.f = m.Approx;
            %[t,y] = m.simulate;
            %m.plot(t,y);
            
            r = m.buildReducedModel;
            ApproxVisualizer(r);
            res = true;
        end
    end
    
end

