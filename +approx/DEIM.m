classdef DEIM < approx.BaseApprox
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
        fOrder = 10;
        
        % The full approximation base
        u;
        
        % Interpolation points
        pts;
        
        % The U matrix for the current Order.
        U;
        
        f;
        
        jrow;
        
        jend;
        
        % The x-component selection matrix
        S;
    end
    
    methods
        function this = DEIM
            this = this@approx.BaseApprox;
            this.CustomProjection = true;
            %this.JSparsityPattern
            this.TimeDependent = true;
            %this.CustomJacobian = true;
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
            
            this.MultiArgumentEvaluations = this.f.MultiArgumentEvaluations;
            
            % Trigger computation of U matrix etc
            this.Order = this.fOrder;
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            fx = this.U * this.f.evaluateComponents(this.pts(1:this.fOrder),...
                this.jend(1:this.fOrder),...
                this.S*x, t, mu);
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected = project@approx.BaseApprox(this, V, W, projected);
            projected.updateOrderData;
        end
        
        function copy = clone(this)
            copy = approx.DEIM;
            copy = clone@approx.BaseApprox(this, copy);
            copy.fOrder = this.fOrder;
            copy.u = this.u;
            copy.pts = this.pts;
            copy.U = this.U;
            copy.f = this.f;
            copy.jrow = this.jrow;
            copy.jend = this.jend;
            copy.S = this.S;
            copy.MaxOrder = this.MaxOrder;
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
            o = this.fOrder;
            jr = [];
            this.jend = zeros(1,o);
            for i=1:this.fOrder
                jr = [jr this.f.getComponentArgumentIndices(this.pts(i))];%#ok
                this.jend(i) = length(jr);
            end
            this.jrow = jr;
            
            n = size(this.u,1);
            P = sparse(this.pts(1:o),1:o,ones(o,1),n,o);
            this.U = this.u(:,1:this.fOrder) * inv(P'*this.u(:,1:this.fOrder));
            if ~isempty(this.W)
                this.U = this.W'*this.U;
            end
            % Compose x-entry selection matrix
            len = this.jend(o);
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
        function m = test_DEIM
            m = models.pcd.PCDModel(1);
            m.EnableTrajectoryCaching = false;
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 40;
            m.System.Params(1).Desired = 10;
            m.SpaceReducer = spacereduction.PODGreedy;
            m.off1_createParamSamples;
            m.off2_genTrainingData;
            m.off3_computeReducedSpace;
            m.off4_genApproximationTrainData;
            save DEIM m;
        end
    end
    
end