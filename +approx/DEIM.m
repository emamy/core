classdef DEIM < approx.BaseApprox
% DEIM: Implements the DEIM-Algorithm from [CS09]
%
% [CS09] Chaturantabut, S. & Sorensen, D.
% Discrete Empirical Interpolation for nonlinear model reduction
% Proc. of CDC/CCC 2009, 2009, pp. 4316 -4321
%
% Additionally allows to estimate the DEIM error a-posteriori by using the 
% next `M'` orders of the current DEIM approximation.
%
% @author Daniel Wirtz @date 2012-03-26
%
% @change{0,6,dw,2012-05-29} Besides many work-in-progress changes, now the
% evaluate function is implemented directly. further, some matrices needed
% for the DEIMEstimator have been included and are being computed.
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% @todo think of exporting the jrow, jend, S properties to ACompEvalCoreFun
% and precompute stuff there; projection happens at that stage
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
        % The actual order `M` for the current DEIM approximation.
        %
        % As second element, the order `M'` of the error estimation DEIM
        % approximation can be specified. If not given, `M'=0` is
        % automatically set.
        %
        % See also: MaxOrder
        %
        % @type rowvec<integer> @default [10 4]
        Order;
    end
    
    properties%(Access=private)
        fOrder = [10 4];
        
        % The full approximation base
        u;
        
        % Interpolation points
        pts;
        
        % The U matrix for the current Order.
        U;
        
        % Some matrices for M+M' error estimation
        Uerr1;
        Uerr2;
        
        M1;
        M2;
        
        % Extra matrix used in the error.DEIMEstimator class
        Uerr1_VW;
        
        % The function which DEIM is applied to
        %
        % Is a subclass of dscomponents.ACompEvalCoreFun
        f;
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
            this.computeDEIM(model.System.f, model.Data.ApproxTrainData.fxi);
        end
        
        function computeDEIM(this, f, fxi)
            this.f = f;
            if ~isa(this.f,'dscomponents.ACompEvalCoreFun');
                error('Cannot use DEIM with non ACompEvalCoreFun-implementing functions.');
            end
            
            %% Generate u_1 ... u_m base
            p = general.POD;
            p.Mode = 'abs';
            p.Value = this.MaxOrder;
            p.UseSVDS = size(fxi,1) > 10000;
            this.u = p.computePOD(fxi);
            
            % Get interpolation points
            this.pts = this.getInterpolationPoints(this.u);
            
            this.MultiArgumentEvaluations = this.f.MultiArgumentEvaluations;
            
            % Trigger computation of U matrix etc
            this.Order = this.fOrder;
        end
        
        function fx = evaluate(this, x, t, mu)
            fx = this.U * this.f.evaluateComponentSet(1, x, t, mu);
        end
        
        function fx = evaluateCoreFun(varargin)
            % do nothing as evaluate is overridden directly
        end
        
        function err = getEstimatedError(this, x, t, mu)
            if this.fOrder(2) == 0
                error('No error estimation possible with zero ErrorOrder property');
            end
            err = this.Uerr1 * this.f.evaluateComponentSet(1, x, t, mu) ...
                - this.Uerr2 * this.f.evaluateComponentSet(2, x, t, mu);
        end
        
        function [res, pm] = computeDEIMErrors(this, atd, orders, errorders)
            oldo = this.fOrder;
            if nargin < 4
                if nargin < 3
                    orders = 1:this.MaxOrder-1;
                end
                errorders = [];
                neworders = [];
                no = length(orders);
                for i=1:no
                    new = 1:(this.MaxOrder-orders(i));
                    errorders = [errorders new];%#ok
                    neworders = [neworders orders(i)*ones(1,length(new))];%#ok
                end
                orders = neworders;
            elseif length(orders) ~= length(errorders)
                error('If both the orders and error orders are given they must have same number of elements');
            end
            no = length(orders);
            % State space error function
            efun = @Norm.L2;
            % Elements error functions
            sumfun = @(x)Norm.Linf(x');
            fxinorm = efun(atd.fxi);
            
            res = zeros(6,no);
%             res = zeros(8,no);
            pi = tools.ProcessIndicator(sprintf('Computing DEIM errors and estimates for %d Order/ErrOrder settings',no),no);
            co = [];
            for i = 1:no
                o = orders(i);
                eo = errorders(i);
                this.Order = [o eo];
                res(1:2,i) = [o; eo];
                
                % Absolute true error (only for new orders)
                if isempty(co) || o ~= co
                    afxi = this.evaluateCoreFun(atd.xi,atd.ti,atd.mui);
                    hlp = efun(atd.fxi-afxi);
                    res(3,i) = sumfun(hlp);
                    % Relative true error 
                    res(4,i) = sumfun(hlp./fxinorm);
                    co = o;
                else
                    res(3:4,i) = res(3:4,i-1);
                end
                % Estimated absolute/rel errors
                hlp = efun(this.getEstimatedError(atd.xi,atd.ti,atd.mui));
                res(5,i) = sumfun(hlp);
                res(6,i) = sumfun(hlp./fxinorm);
                
%                 % Compute actual error between M and M' approximations
%                 this.Order = [sum(this.fOrder) 0];
%                 % Get order+errorder eval
%                 hlp = efun(this.evaluateCoreFun(atd.xi,atd.ti,atd.mui) - afxi);
%                 res(7,i) = sumfun(hlp);
%                 res(8,i) = sumfun(hlp./fxinorm);
                pi.step;
            end
            pi.stop;
            if nargout == 2
                pm = approx.DEIM.plotDEIMErrs(res);
            end
            this.Order = oldo;
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected = project@approx.BaseApprox(this, V, W, projected);
            % Important: Project the component evaluable function, too!
            projected.f = this.f.project(V,W);
            projected.updateOrderData;
        end
        
        function copy = clone(this)
            copy = approx.DEIM;
            copy = clone@approx.BaseApprox(this, copy);
            copy.fOrder = this.fOrder;
            copy.u = this.u;
            copy.pts = this.pts;
            copy.U = this.U;
            copy.f = this.f.clone;
            copy.MaxOrder = this.MaxOrder;
            copy.Uerr1 = this.Uerr1;
            copy.Uerr2 = this.Uerr2;
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
                fprintf('DEIM interpolation points [%s] with values [%s]\n',...
                    general.Utils.implode(pts,' ','%d'),general.Utils.implode(v,' ','%2.2e'));
            end
        end
        
        function updateOrderData(this)
            o = this.fOrder(1);
            n = size(this.u,1);
            
            P = sparse(this.pts(1:o),1:o,ones(o,1),n,o);
            this.U = this.u(:,1:o) / (P'*this.u(:,1:o));
            
            % Set primary point set in ACompEvalCoreFun
            this.f.setPointSet(1,this.pts(1:o));
            
            om = this.fOrder(2);
            if om > 0
                % Use second point set in ACompEvalCoreFun
                this.f.setPointSet(2,this.pts(o+1:o+om));
                
                Perr = sparse(this.pts(o+1:o+om),1:om,ones(om,1),n,om);
                
                Um = this.u(:,1:o);
                Umd = this.u(:,o+1:o+om);
                A = P'*Um;
                B = P'*Umd;
                C = Perr'*Um;
                D = Perr'*Umd;
                E = C / A;
                F = D - E*B;
                this.Uerr2 = ((Um/A)*B - Umd) / F;
                this.Uerr1 = this.Uerr2 * E;
            end
            
            if ~isempty(this.W)
                this.U = this.W'*this.U;
                if om > 0
                    % Compute values for error estimator
                    this.M1 = this.Uerr1 + (Um - this.V*(this.W'*Um))/A;
                    this.M2 = this.Uerr2;
                    
                    % Project DEIM error estimation values
                    this.Uerr1 = this.W'*this.Uerr1;
                    this.Uerr2 = this.W'*this.Uerr2;
                end
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
            if numel(value) ~= 2
                if numel(value) ~= 1
                    error('Order must be either a scalar or two element integer vector');
                else
                    value(2) = 0;
                end
            end
            if any(value) > size(this.u,2) || any(value) < 1
                error('Invalid Order/ErrOrder value. Allowed are integers in [1, %d]',size(this.u,2));
            elseif sum(value) > this.MaxOrder
                error('Order (%d) and ErrOrder (%d) values may not exceed MaxOrder (%d)',value,this.MaxOrder);
            end
            this.fOrder = value;
            
            this.updateOrderData;
        end
    end
    
    methods(Static)
%         function [res, pm] = computeOrderApproxOnATD(model)
%             a = model.Approx;
%             atd = model.Data.ApproxTrainData;
%             res = zeros(2,a.MaxOrder);
%             for i=1:a.MaxOrder
%                 a.Order = i;
%                 afxi = a.evaluate(atd.xi,atd.ti,atd.mui);
%                 res(1,i) = mean(Norm.L2(atd.fxi-afxi));
%                 res(2,i) = mean(Norm.L2(atd.fxi-afxi)./Norm.L2(atd.fxi));
%             end
%             pm = tools.PlotManager(false,1,2);
%             pm.nextPlot('mean_abs','Average absolute approximation error over all DEIM-orders on ApproxTrainData','order','mean abs error on training data');
%             semilogy(res(1,:));
%             pm.nextPlot('mean_rel','Average relative approximation error over all DEIM-orders on ApproxTrainData','order','mean rel error on training data');
%             semilogy(res(2,:));
%             pm.done;
%         end
        
        function pm = plotDEIMErrs(res, pm)
            %% Plotting
            tri = delaunay(res(1,:),res(2,:));
            if nargin < 2
                pm = tools.PlotManager(false,2,3);
                pm.FilePrefix = 'deim';
            end
            h = pm.nextPlot('true_abs','Linf-L2 absolute error','order','error order');
            doplot(h,tri,res(1,:),res(2,:),res(3,:));
            h = pm.nextPlot('true_rel','Linf-L2 relative  error','order','error order');
            doplot(h,tri,res(1,:),res(2,:),res(4,:));
            
            h = pm.nextPlot('est_abs','Estimated absolute error','order','error order');
            doplot(h,tri,res(1,:),res(2,:),res(5,:));
            h = pm.nextPlot('est_rel','Estimated relative error','order','error order');
            doplot(h,tri,res(1,:),res(2,:),res(6,:));
                        
            h = pm.nextPlot('abs_diff','|true - estimated| absolute error','order','error order');
            doplot(h,tri,res(1,:),res(2,:),abs(res(5,:)-res(3,:)));
            view(0,90);
            h = pm.nextPlot('rel_diff','|true - estimated| relative error','order','error order');
            doplot(h,tri,res(1,:),res(2,:),abs(res(6,:)-res(4,:)));
            view(71,52);
            
%             h = pm.nextPlot('abs_diff','|true - dir. est.| absolute error','order','error order');
%             doplot(h,tri,res(1,:),res(2,:),abs(res(7,:)-res(3,:)));
%             h = pm.nextPlot('rel_diff','|true - dir. est.| relative error','order','error order');
%             doplot(h,tri,res(1,:),res(2,:),abs(res(8,:)-res(4,:)));
%             
%             h = pm.nextPlot('abs_diff','|dir est - ref. est.| absolute error','order','error order');
%             doplot(h,tri,res(1,:),res(2,:),abs(res(5,:)-res(7,:)));
%             h = pm.nextPlot('rel_diff','|dir est - ref. est.| relative error','order','error order');
%             doplot(h,tri,res(1,:),res(2,:),abs(res(6,:)-res(8,:)));
%             
%             h = pm.nextPlot('est_abs','Direct estimation of absolute error','order','error order');
%             doplot(h,tri,res(1,:),res(2,:),res(7,:));
%             h = pm.nextPlot('est_rel','Direct estimation of relative error','order','error order');
%             doplot(h,tri,res(1,:),res(2,:),res(8,:));
            
            pm.done;
            
            function p = doplot(h, tri, x, y, z)
                p = tools.LogPlot.logtrisurf(h, tri, x, y, z);
                hold on;
                [~, hc] = tricontour([x; y]', tri, z', get(h,'ZTick'));
                set(hc,'EdgeColor','k');
                hold off;
                view(0,0);
            end
        end
        
        function [m, r] = test_DEIM
            m = models.pcd.PCDModel(1);
            m.EnableTrajectoryCaching = false;
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 40;
            m.System.Params(1).Desired = 10;
            m.SpaceReducer = spacereduction.PODGreedy;
%             m.offlineGenerations;
            m.off1_createParamSamples;
            m.off2_genTrainingData;
            m.off3_computeReducedSpace;
            m.off4_genApproximationTrainData;
            m.off5_computeApproximation;
            m.Approx.Order = [20 4];
            r = m.buildReducedModel;
            save DEIM m r;
        end
    end
    
end