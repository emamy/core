classdef KernelEI < approx.BaseApprox
% KernelEI: DEIM approximation using kernel expansions for
% function/operator evaluations
%
% Implements the ideas described in @cite W13 .
%
% @author Daniel Wirtz @date 2012-03-26
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetObservable)
        % The maximum order up to which the DEIM approximation should be computed.
        %
        % This corresponds to the maximum number `m` that can be choosen as approximation order.
        %
        % @type integer @default 40
        MaxOrder = 40;
        
        % The variant to use.
        %
        % Possible values are 1 and 2. 1 means that all components are
        % interpolated simultaneously and only the current KernelEI.Order
        Variant = 2;
        
        % An approx.ABase approximation algorithm that is used to learn the
        % component functions (either simultaneously or component-wise
        % depending on Variant)
        %
        % @type approx.ABase @default []
        %
        % See also: Variant
        Algorithm;
    end
    
    properties(SetAccess=private)
        % The Kernel expansion computed if using Variant 1.
        %
        % @type kernels.KernelExpansion
        V1Expansion = [];
        
        % The array of kernel expansions of using Variant 2.
        %
        % @type rowvec<kernels.KernelExpansion>
        V2Expansions = {};
    end

    properties(SetObservable, Dependent)
        % The actual order for the current DEIM approximation.
        %
        % See also: MaxOrder
        %
        % @type integer @default 10
        Order;
    end

    properties(Access=private)
        fOrder = 5;
        
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
        end
        
        function approximateSystemFunction(this, model)
            this.f = model.System.f;
            if ~isa(this.f,'dscomponents.ACompEvalCoreFun');
                error('Cannot use DEIM with no ACompEvalCoreFun core functions.');
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
            SP = this.f.JSparsityPattern;
            jr = [];
            je = zeros(1,length(this.pts));
            for i=1:length(this.pts)
                sprow = SP(this.pts(i),:);
                inew = find(sprow);
                jr = [jr inew];%#ok
                je(i) = length(jr);
            end
            this.jrow = jr;
            this.jend = je;
            
            % Trigger computation of U matrix etc
            this.Order = this.fOrder;
            
            this.V1Expansion = [];
            this.V2Expansions = [];
            if this.Variant == 1
                % Compute Vz and f(Vz) values globally
                zi = model.Data.W'*atd.xi;
                fzi = model.System.f.evaluate(model.Data.V*zi,atd.ti,atd.mui);
                vatd = data.ApproxTrainData(zi,atd.ti,atd.mui);
                mo = this.MaxOrder;
                vatd.fxi = sparse(this.pts,1:mo,ones(mo,1),size(atd.fxi,1),mo)'*fzi;
                this.V1Expansion = this.Algorithm.computeApproximation(vatd);
            elseif this.Variant == 2
                % Compose x-entry selection matrix
                pi = ProcessIndicator(sprintf('KernelEI: Computing %d approximations',this.MaxOrder),this.MaxOrder);
                for idx = 1:this.MaxOrder
                    % Select the elements of x that are effectively used in f
                    off = 0;
                    if idx > 1
                        off = this.jend(idx-1);
                    end
                    xireq = atd.xi(this.jrow(off+1:this.jend(idx)),:);
                    latd = data.ApproxTrainData(xireq, atd.ti, atd.mui);
                    latd.fxi = atd.fxi(this.pts(idx),:);
                    this.V2Expansions{idx} = this.Algorithm.computeApproximation(latd);
                    pi.step;
                end
                pi.stop;
            end
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            if this.Variant == 1
                fu = this.V1Expansion.evaluate(x, t, mu);
                c = fu(1:this.fOrder);
            elseif this.Variant == 2
                c = zeros(this.fOrder,size(x,2));
                for o = 1:this.fOrder
                    % Select the elements of x that are effectively used in f
                    % S is already of the size of all required elements, so to this.jrow
                    % indexing is required
                    off = 1;
                    if o > 1
                        off = off + this.jend(o-1);
                    end
                    c(o,:) = this.V2Expansions{o}.evaluate(this.S(off:this.jend(o),:)*x,t,mu);
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
            if ~isempty(this.V1Expansion)
                copy.V1Expansion = this.V1Expansion.clone;
            end
            copy.Variant = this.Variant;
            copy.jend = this.jend;
            copy.jrow = this.jrow;
            n = length(this.V2Expansions);
            for i=1:n
                copy.V2Expansions{i} = this.V2Expansions{i}.clone;
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
                    Utils.implode(pts,' ','%d'),Utils.implode(v,' ','%2.2e'));
            end
        end
        
        function updateOrderData(this)
            n = size(this.u,1);
            o = this.fOrder;
            P = sparse(this.pts(1:o),1:o,ones(o,1),n,o);
            this.U = this.u(:,1:this.fOrder) * ...
                inv(P'*this.u(:,1:this.fOrder));
            if ~isempty(this.W)
                this.U = this.W'*this.U;
            end
            
            len = this.jend(end);
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
        function res = test_KernelEI
            m = models.pcd.PCDModel(1);
            m.dt = .1;
            m.T = 1;
            m.EnableTrajectoryCaching = false;
            
            % Define prototype expansion
            ec = kernels.config.ParamTimeExpansionConfig;
            ec.StateConfig = kernels.config.GaussConfig('G',1:3);
            ec.TimeConfig = [];
            ec.ParamConfig = kernels.config.GaussConfig('G',1:3);
            
            a = approx.algorithms.VKOGA;
            a.MaxExpansionSize = 300;
            a.UsefScaling = false;
            a.UsefPGreedy = false;
            a.ExpConfig = ec;
            
            kei = approx.KernelEI;
            kei.Algorithm = a;
            kei.Variant = 2;
            
            m.Approx = kei;
            m.Approx.MaxOrder = 5;
            m.System.Params(1).Desired = 2;
            m.SpaceReducer = spacereduction.PODGreedy;
            m.SpaceReducer.Eps = 1e-2;
            m.offlineGenerations;
            
            mu = m.getRandomParam;
            r = m.buildReducedModel;
            r.simulate(mu);
            
            kei.Variant = 1;
            m.off5_computeApproximation;
            
            r = m.buildReducedModel;
            r.simulate(mu);

            res = true;
        end
    end
    
end