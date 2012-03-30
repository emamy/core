classdef KernelEI < approx.BaseApprox
% DEIM: 
%
%
%
% @author Daniel Wirtz @date 2012-03-26
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The maximum order up to which the DEIM approximation should be computed.
        %
        % This corresponds to the maximum number `m` that can be choosen as approximation order.
        %
        % @type integer @default 40
        MaxOrder = 40;
        
        kexp;
        
        exps;
        
        variant = 2;
    end

    properties(SetObservable, Dependent)
        % The actual order for the current DEIM approximation.
        %
        % See also: MaxOrder
        %
        % @type integer @default 10
        Order;
    end
%     
    properties%(Access=private)
        fOrder = 10;
        
        % The full approximation base
        u;
        
        % Interpolation points
        pts;
        
        % The U matrix for the current Order.
        U;
        
        S;
        
        f;
        
        jrow;
        
        jend;
    end
    
    methods
        function this = KernelEI
            this = this@approx.BaseApprox;
            this.CustomProjection = true;
            this.MultiArgumentEvaluations = true;
            this.TimeDependent = true;
            %this.CustomJacobian = true;
            %this.JSparsityPattern
        end
        
        function approximateSystemFunction(this, model)
            this.f = model.System.f;
            if ~isa(this.f,'dscomponents.IComponentEvaluable');
                error('Cannot use DEIM with no IComponentEvaluable core functions.');
            end
            
            atd = model.Data.ApproxTrainData;
            
            %% Generate u_1 ... u_m base
            p = general.POD;
            p.Mode = 'abs';
            p.Value = this.MaxOrder;
            p.UseSVDS = size(atd.fxi,1) > 10000;
            this.u = p.computePOD(atd.fxi);
            
            % Get interpolation points
            this.pts = this.getInterpolationPoints(this.u);
            
            % Compose argument indices arrays
            jr = [];
            this.jend = zeros(1,this.MaxOrder+1);
            this.jend(1) = 0;
            for i=1:this.MaxOrder
                jr = [jr this.f.getComponentArgumentIndices(this.pts(i))];%#ok
                this.jend(i+1) = length(jr);
            end
            this.jrow = jr;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.NumGammas = 5;
            %a.Dists = 2;
            a.MaxExpansionSize = 100;
            a.UsefScaling = false;
            a.UseOGA = false;
            
            k = kernels.ParamTimeKernelExpansion;
            k.Kernel = kernels.GaussKernel;
            if ~isempty(atd.ti)
                k.TimeKernel = kernels.GaussKernel;
            end
            if ~isempty(atd.mui)
                k.ParamKernel = kernels.GaussKernel;
            end
            
            % Trigger computation of U matrix etc
            this.Order = this.fOrder;
            
            if this.variant == 1
                % Compute Vz and f(Vz) values
                zi = model.Data.W'*atd.xi;
                fzi = model.System.f.evaluate(model.Data.V*zi,atd.ti,atd.mui);
                vatd = data.ApproxTrainData(zi,atd.ti,atd.mui);
                mo = this.MaxOrder;
                vatd.fxi = sparse(this.pts,1:mo,ones(mo,1),size(atd.fxi,1),mo)'*fzi;
                this.kexp = k;
                a.computeApproximation(this.kexp, vatd);
            elseif this.variant == 2
                % Compose x-entry selection matrix
                this.exps = kernels.ParamTimeKernelExpansion.empty(0,this.MaxOrder);
                pi = tools.ProcessIndicator(sprintf('KernelEI: Computing %d approximations',this.MaxOrder),this.MaxOrder);
                for idx = 1:this.MaxOrder
                    % Select the elements of x that are effectively used in f
                    xireq = atd.xi(this.jrow(this.jend(idx)+1:this.jend(idx+1)),:);
                    latd = data.ApproxTrainData(xireq, atd.ti, atd.mui);
                    latd.fxi = atd.fxi(this.pts(idx),:);
                    a.computeApproximation(k, latd);
                    this.exps{idx} = k.clone;
                    pi.step;
                end
                pi.stop;
            end
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            if this.variant == 1
                fu = this.kexp.evaluate(x, t, mu);
                c = fu(1:this.fOrder);
            elseif this.variant == 2
                c = zeros(this.fOrder,size(x,2));
                for o = 1:this.fOrder
                    % Select the elements of x that are effectively used in f
                    % S is already of the size of all required elements, so to this.jrow
                    % indexing is required
                    xidx = (this.jend(o)+1):this.jend(o+1);
                    c(o,:) = this.exps{o}.evaluate(this.S(xidx,:)*x,t,mu);
                end
            end
            fx = this.U * c;
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected = project@approx.BaseApprox(this, V, W, projected);
            projected.updateOrderData;
        end
        
        function copy = clone(this)
            copy = approx.KernelEI;
            copy = clone@approx.BaseApprox(this, copy);
            copy.fOrder = this.fOrder;
            copy.u = this.u;
            copy.pts = this.pts;
            copy.U = this.U;
            copy.f = this.f;
            copy.jrow = this.jrow;
            copy.jend = this.jend;
            copy.MaxOrder = this.MaxOrder;
            if ~isempty(this.kexp)
                copy.kexp = this.kexp.clone;
            end
            copy.variant = this.variant;
            copy.jend = this.jend;
            copy.jrow = this.jrow;
            n = length(this.exps);
            copy.exps = kernels.ParamTimeKernelExpansion.empty(0,n);
            for i=1:n
                copy.exps{i} = this.exps{i}.clone;
            end
            copy.S = this.S;
        end
    end
    
    methods(Access=private)
        function pts = getInterpolationPoints(~, u)
            n =size(u,1);
            m = size(u,2);
            pts = zeros(1, m);
            v = pts;
            [v(1), pts(1)] = max(abs(u(:,1)));
            P = zeros(n,1);
            P(pts(1)) = 1;
            for i=2:m
                c = (P'*u(:,1:(i-1))) \ (P'*u(:,i));
                [v(i), pts(i)] = max(abs(u(:,i) - u(:,1:(i-1))*c));
                P = sparse(pts(1:i),1:i,ones(i,1),n,i);
            end
            if KerMor.App.Verbose > 2
                fprintf('DEIM interpolation points [%s] with values [%s]',...
                    general.Utils.implode(pts,' ','%d'),general.Utils.implode(v,' ','%2.2e'));
            end
        end
        
        function updateOrderData(this)
            n = size(this.u,1);
            o = this.fOrder;
            P = sparse(this.pts(1:o),1:o,ones(o,1),n,o);
            this.U = this.u(:,1:this.fOrder) * inv(P'*this.u(:,1:this.fOrder));
            if ~isempty(this.W)
                this.U = this.W'*this.U;
            end
            
            len = this.jend(this.fOrder+1);
            sel = this.jrow(1:len);
            if ~isempty(this.V)
                this.S = this.V(sel,:);
            else
                this.S = sparse(1:len,sel,ones(len,1),len,n);
            end
        end
    end
    
    %% Getter & Setter
    methods
        function o = get.Order(this)
            o = this.fOrder;
        end
        
        function set.Order(this, value)
            if isempty(this.u) || isempty(this.pts)
                error('Cannot set DEIM order as approximateSystemFunction has not been called yet.');
            end
            if value > size(this.u,2) || value < 1
                error('Invalid Order value. Allowed are integers in [1, %d]',size(this.u,2));
            end
            this.fOrder = value;
            
            this.updateOrderData;
        end
    end
    
    methods(Static)
        function m = test_KernelEI
            m = models.pcd.PCDModel(1);
            m.EnableTrajectoryCaching = false;
            m.Approx = approx.KernelEI;
            m.Approx.MaxOrder = 30;
            m.System.Params(1).Desired = 10;
            m.SpaceReducer = spacereduction.PODGreedy;
            m.off1_createParamSamples;
            m.off2_genTrainingData;
            m.off3_computeReducedSpace;
            m.off4_genApproximationTrainData;
            save KernelEI m;
        end
    end
    
end