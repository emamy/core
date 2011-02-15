classdef TPWLApprox < approx.BaseApprox
    %TPWLAPPROX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The weight function `w`
        WeightFun;
        
        % The distance within there has to be an expansion point for each
        % x. Gets multiplied by `\sqrt{d}`, where `d` denotes the spatial
        % dimension of the snapshots (projection training data)
        EpsRad = .05;
        
        % The minimum value a weight needs to have after normalization in
        % order to affect the evaluation. Set this value to a certain
        % threshold in order to restrict globally supported weight
        % functions (i.e. the Gaussian) to a local support.
        MinWeightValue = 1e-10;
    end
    
    properties(Access=private)
        % The gradient matrix
        Ai;
        bi;
        xi;
    end
    
    methods
        
        function this = TPWLApprox
            this.WeightFun = kernels.GaussKernel(25);
        end
        
        function copy = clone(this)
            % Implements ICloneable.clone()
            copy = approx.TPWLApprox;
            copy = clone@approx.BaseApprox(this, copy);
            copy.Ai = this.Ai;
            copy.bi = this.bi;
            copy.xi = this.xi;
            copy.EpsRad = this.EpsRad;
            copy.MinWeightValue = this.MinWeightValue;
            if isa(this.WeightFun,'ICloneable')
                copy.WeightFun = this.WeightFun.clone;
            else
                copy.WeightFun = this.WeightFun;
            end
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
            xsel = reshape(meshgrid(1:as,1:size(x,2)),1,[]);
            diff = repmat(x,1,as)-this.xi(:,xsel);
            w = this.WeightFun.evaluateScalar(sqrt(sum(diff.^2,1)));
            
            y = zeros(size(x));
            for i=1:as
                % Extract the weights for x and the current xi
                wl = w((i-1)*size(x,2)+1:i*size(x,2));
                % Normalize local weights
                wl = wl / sum(wl);
                % Find relevant weights
                rel = find(wl > this.MinWeightValue);
                % Only update entries for relevant weights
                y(:,rel) = y(:,rel) ...
                    +( this.Ai{i}*x(:,rel)...
                       +this.bi(:,i*ones(1,length(rel)))...
                     )*diag(wl(rel));
            end
        end
        
        function atd = selectTrainingData(this, modeldata)
            sn = modeldata.ProjTrainData;
            x = sn(4:end,:);
            d = sqrt(size(x,1))*this.EpsRad;
            selidx = 1;
            idx = 1; 
            cur = x(:,idx);
            while idx < size(x,2)
                while(idx < size(x,2) && norm(x(:,idx)-cur) < d)
                    idx = idx+1;
                end
                selidx(end+1) = idx;%#ok
                cur = x(:,idx);
                idx = idx+1;
            end
            atd = sn(:,selidx);
        end
    end
    
    methods(Access=protected)
        function gen_approximation_data(this, model, xi, ti, mui)
            % Implements BaseApprox abstract template method
            %
            % @todo parallelize
            
            as = size(xi,2);
            N = size(xi,1);
            h = 0.01;
            
            this.Ai = cell(0,as);
            dh = diag(ones(1,N)*h);
            for i = 1:as
                sel = ones(1,N)*i;
                xipt = xi(:,sel);
                if isempty(mui); mu = []; else mu = mui(:,sel); end
                tmp = (model.Data.ApproxfValues(:,sel)-...
                    model.f.evaluate(xipt+dh,ti(sel),mu))/h;
                this.Ai{i} = tmp;
                this.bi(:,i) = model.Data.ApproxfValues(:,i) - tmp*xi(:,i);
            end
            this.xi = xi;
            
%             this.Ai = sparse(N*as,N*as);%xi(:,idx) + repmat(diag(ones(1,N))*h,1,as);
%             dh = diag(ones(1,N)*h);
%             for i = 1:as
%                 % compute jacobian at xi
%                 %fprintf('Iteration %d\n',i);
%                 %for j=1:N
%                 sel = ones(1,N)*i;
%                 xipt = xi(:,sel);
%                 fxipt = model.Data.ApproxfValues(:,sel);
%                 if isempty(mui); mu = []; else mu = mui(:,sel); end
%                 tmp = (fxipt-model.f.evaluate(xipt+dh,ti(sel),mu))/h;
%                 this.Ai((i-1)*N+1:i*N,(i-1)*N+1:i*N) = tmp;
%             end
%             this.xi = xi;
%             % store as large column vectors
%             this.bi = reshape(model.Data.ApproxfValues,[],1) - this.Ai*reshape(xi,[],1);
        end
    end
    
    methods(Static)
        function res = test_TWPLApprox
            dims = 100;
            epsrad = .05;
            
            m = models.synth.KernelTest(dims,false);
            m.Approx = approx.TPWLApprox;
            m.Approx.EpsRad = epsrad;
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
            
            r.analyze;
        end
    end
    
end

