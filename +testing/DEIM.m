classdef DEIM
% DEIM: Tests regarding the DEIM method.
%
% See also: approx.DEIM
%
% @author Daniel Wirtz @date 2012-05-03
%
% @new{0,6,dw,2012-05-03} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        function analysis_DEIM_approx(m)
            ma = tools.ModelAnalyzer(m.buildReducedModel);
            
            mu = m.Data.ParamSamples(:,5);
%             [t, pm] = ma.compareRedFull(mu);
%             t.Format = 'tex';
%             t.saveToFile('comp_red_full.tex');
%             pm.savePlots('.',{'fig','jpg','eps'}, true);
            
            d = m.Approx;
            res = testing.DEIM.computeDEIMErrors(d, m.Data.ApproxTrainData);
            %pm = tools.PlotManager(true);
            pm = tools.PlotManager(false,1,2);
            pm.FilePrefix = 'DEIM_errors';
            testing.DEIM.plotDEIMErrs(res, pm);
            %pm.savePlots('.',{'fig','jpg','eps'}, true);
            pm.savePlots('.',{'fig','jpg'}, true);
            
            [etrue, EE, ED] = ma.getTrajApproxErrorDEIMEstimates(mu,[]);
%             pm = tools.PlotManager(true);
            pm = tools.PlotManager(false,1,2);
            pm.FilePrefix = 'traj_approx_err';
            ma.getTrajApproxErrorDEIMEstimates_plots(m.scaledTimes, etrue, EE, ED, pm);
%             pm.savePlots('.',{'fig','jpg','eps'}, true);
            pm.savePlots('.',{'fig','jpg'}, true);
            
            [minreq, relerrs, orders, t] = ma.getMeanRequiredErrorOrders;
            t.Format = 'tex';
            t.saveToFile('meanrequirederrororders.tex');
            
            save analysis_DEIM_approx;
        end
        
        function [res, pm] = computeDEIMErrors(deim, atd, orders, errorders)
            % Computes the DEIM approximation error over the given training
            % data for specified DEIM orders (M) and DEIM error orders
            % (M').
            %
            % The error over any training data point is measured L2 in
            % state and Linf over all training points.
            oldo = deim.Order;
            if nargin < 4
                if nargin < 3
                    orders = 1:deim.MaxOrder-1;
                end
                errorders = [];
                neworders = [];
                no = length(orders);
                for i=1:no
                    new = 1:(deim.MaxOrder-orders(i));
                    errorders = [errorders new];
                    neworders = [neworders orders(i)*ones(1,length(new))];
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
            
            if isa(deim,'general.MatrixDEIM')
                oldst = deim.Qk;
                deim.setSimilarityTransform([]);
            end
            
            res = zeros(6,no);
%             res = zeros(8,no);
            pi = tools.ProcessIndicator(sprintf('Computing DEIM errors and estimates for %d Order/ErrOrder settings',no),no);
            co = [];
            xi = atd.xi;
            if ~isempty(deim.V)
                xi = deim.V'*xi;
            end
            for i = 1:no
                o = orders(i);
                eo = errorders(i);
                deim.Order = [o eo];
                res(1:2,i) = [o; eo];
                
                % Absolute true error (only for new orders)
                if isempty(co) || o ~= co
                    if isa(deim,'general.MatrixDEIM')
                        for j=1:size(xi,2)
                            afx = deim.evaluate(xi(:,j),atd.ti(j),atd.mui(:,j));
                            hlp(j) = efun(atd.fxi(:,j)-afx(:));
                        end
                    else
                        afxi = deim.evaluate(xi,atd.ti,atd.mui);
                        hlp = efun(atd.fxi-afxi);
                    end
                    
                    res(3,i) = sumfun(hlp);
                    % Relative true error 
                    res(4,i) = sumfun(hlp./fxinorm);
                    co = o;
                else
                    res(3:4,i) = res(3:4,i-1);
                end
                % Estimated absolute/rel errors
                hlp = efun(deim.getEstimatedError(xi,atd.ti,atd.mui));
                res(5,i) = sumfun(hlp);
                res(6,i) = sumfun(hlp./fxinorm);
                
                pi.step;
            end
            pi.stop;
            if nargout == 2
                pm = testing.DEIM.plotDEIMErrs(res);
            end
            deim.Order = oldo;
            if isa(deim,'general.MatrixDEIM')
                deim.setSimilarityTransform(oldst);
            end
        end
        
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
            hold on;
            sh = doplot(h,tri,res(1,:),res(2,:),1e-2*ones(size(res(6,:))),...
                'EdgeColor','none','FaceColor','k');
            alpha(sh,.2);
            hold off;
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
            
            function p = doplot(h, tri, x, y, z, varargin)
                p = tools.LogPlot.logtrisurf(h, tri, x, y, z, varargin{:});
%                 hold on;
%                 [~, hc] = tricontour([x; y]', tri, z', get(h,'ZTick'));
%                 set(hc,'EdgeColor','k');
%                 hold off;
                view(0,0);
            end
        end
               
        function [t, nof, nor, o, pm] = jacobian_tests(m, pm)
            % Tis function computes the jacobians of both the full and DEIM
            % approximated system functions.
            %
            % The maximum relative error for each DEIM order and
            % jacobian at a certain snapshot vector are computed and
            % displayed.
            atd = m.Data.JacobianTrainData;
            n = min(size(atd.xi,2),500);
            sp = logical(m.System.f.JSparsityPattern);
            M = 3;
            o = round(linspace(1,m.Approx.MaxOrder,M));
            t = zeros(M,n);
            nof = t;
            nor = t;
            r = m.buildReducedModel;
            for k = 1:M
                m.Approx.Order = o(k);
                r.System.f.Order = o(k);
                for i = 1:n
                    x = atd.xi(:,i);
                    z = m.Data.W'*x;
                    x = m.Data.V*z;
                    aj = m.Approx.getStateJacobian(x,atd.ti(i),atd.mui(:,i));
                    j = m.System.f.getStateJacobian(x,atd.ti(i),atd.mui(:,i));
                    rj = r.System.f.getStateJacobian(z,atd.ti(i),atd.mui(:,i));
                    hlp = abs((aj-j)./j);
                    hlp(~sp) = 0;
                    t(k,i) = max(hlp(:));
                    nof(k,i) = norm(full(j));
                    nor(k,i) = norm(rj);
                end
            end
            if nargin < 2
                pm = tools.PlotManager(false,3,1);
            end
            h = pm.nextPlot('max_rel_err',sprintf('Maximum relative errors of DEIM jacobian vs. original one\nmasked to sparsity pattern'));
            plot(h,t')
            lbl = arrayfun(@(e)sprintf('%d',e),o,'Uniform',false);
            legend(h, lbl{:});
            
            h = pm.nextPlot('jfull_norm','Full Jacobian norms','snapshot','L2 matrix norm');
            plot(h,nof');
            legend(h, lbl{:});
            
            h = pm.nextPlot('jred_norm','Reduced Jacobian norms','snapshot','L2 matrix norm');
            plot(h,nor');
            legend(h, lbl{:});
           
            pm.done;
        end
        
        function estimator_tests(m)
            r = m.buildReducedModel;
            e = tools.EstimatorAnalyzer(r);
            pm = tools.PlotManager;
            e.start(.01,[],pm);
            pm.done;
%             pm.FilePrefix = class(m);
%             a = KerMor.App;
%             pm.savePlots(fullfile(a.DataStoreDirectory,...
%                 'deim_estimator_tests'),{'fig','jpg'});
        end
        
        function [etrue, EE, ED, pm] = getTrajApproxErrorDEIMEstimates(r, mu, inputidx)
            % [etrue, EE, ED, pm] = getTrajApproxErrorDEIMEstimates(this, mu, inputidx)
            % Computes the DEIM approximation errors for given mu and input.
            %
            % Depending on the current Order `o` of the DEIM approximation,
            % all possible ErrorOrders `eo = o+1 ... N` to the
            % approx.DEIM.MaxOrder `N` are used estimations for the given
            % trajectory computed.
            %
            % Return values:
            % etrue: The true approximation error using the full
            % trajectory.
            % EE: The matrix of estimated errors, containing the estimated
            % state-space L2 error (columns) for each possible ErrorOrder
            % (rows).
            % time step.
            % ED: The absolute difference between true and estimated error
            % for each ErrorOrder (rows) and timestep (columns)
            % pm: The PlotManager used to create the result plots. If not
            % requested as output, no plotting will be done. @optional
            fm = r.FullModel;
            if nargin < 3
                inputidx = fm.TrainingInputs;
                if nargin < 2
                    mu = fm.Data.ParamSamples;
                end
            end
            f = fm.Approx;
            olde = f.Order(2); % store old setting
            o = f.Order(1);
            mo = f.MaxOrder;
            
            num_mu = max(1,size(mu,2));
            num_in = max(1,length(inputidx));
            etrue = zeros(num_mu*num_in,length(fm.Times));
            EE = zeros(mo-o,length(fm.Times),num_mu*num_in);
            ED = EE;
            
            ma = tools.ModelAnalyzer(r);
            
            % Assume no parameters or inputs
            curmu = [];
            curin = [];
            for eo = 1:mo-o
                f.Order = [o eo];
                % Iterate through all input functions
                cnt = 1;
                for inidx = 1:num_in
                    % Iterate through all parameter samples
                    for muidx = 1:num_mu
                        % Check for inputs
                        if ~isempty(inputidx)
                            curin = inputidx(inidx);
                        end
                        % Check for parameters
                        if size(mu,2) > 0
                            curmu = mu(:,muidx);
                        end

                        [curetrue, ~, t, x] = ma.getTrajApproxError(curmu, curin);
                        eest = Norm.L2(f.getEstimatedError(x, t, repmat(curmu,1,length(t))));
                        EE(eo,:,cnt) = eest;
                        ED(eo,:,cnt) = abs(eest - curetrue);
                        etrue(cnt,:) = curetrue;
                        cnt = cnt+1;
                    end
                end
            end
            f.Order = [o olde]; % restore old setting
            
            if nargout == 4
                pm = this.getTrajApproxErrorDEIMEstimates_plots(r, etrue(1,:), EE(:,:,1), ED(:,:,1));
            end
        end
        
        function pm = getTrajApproxErrorDEIMEstimates_plots(r, etrue, EE, ED, pm)
            % Visualizes the results of the
            % tools.ModelAnalyzer.getTrajApproxErrorDEIMEstimates method.
            f = r.FullModel.Approx;
            o = f.Order(1);
            mo = f.MaxOrder;
            if nargin < 6
                pm = tools.PlotManager(false,2,2);
            end
            t = this.rm.Times;
            h = pm.nextPlot('true_err','True approximation error on trajectory','time','error');
            semilogy(h,t,etrue);
            h = pm.nextPlot('est_err',sprintf('Estimated approximation errors on trajectory\nDEIM order = %d / max %d',f.Order(1),f.MaxOrder),...
                'time', 'DEIM error order');
            [T,O] = meshgrid(t,1:mo-o);
            tools.LogPlot.logsurf(h,T,O,EE,'EdgeColor','none');
            view(-45,45);
            
            h = pm.nextPlot('abs_diff',sprintf('Difference true to estimated error\nDEIM order = %d / max %d',f.Order(1),f.MaxOrder),...
                'time','DEIM error order');
            tools.LogPlot.logsurf(h,T,O,ED,'EdgeColor','none');
            view(-135,45);
            
            h = pm.nextPlot('rel_diff',sprintf('Relative error between true to estimated error\nDEIM order = %d / max %d',f.Order(1),f.MaxOrder),...
                'time','DEIM error order');
            tools.LogPlot.logsurf(h,T,O,ED./repmat(etrue,mo-o,1),'EdgeColor','none','FaceColor','interp');
            hold on;
            p = surf(h,T,O,-ones(size(T))*2,'EdgeColor','red','FaceColor','k');
            alpha(p,.2);
            view(-135,45);
            hold off;
            pm.done;
        end
        
        function [res, relerrs, orders, t] = getMeanRequiredErrorOrders(r, relerrs, orders, mu, inidx)
            % [res, relerrs, orders, t] = getMeanRequiredErrorOrders(this, relerrs, orders, mu)
            % This high-level function computes the minimum values of M' in
            % order to have the relative error of true to estimated DEIM
            % error smaller than the specified relerrs.
            %
            % Parameters:
            % relerrs: The relative errors that may occur between true and
            % estimated approximation error, maximized over time of each
            % trajectory, averaged over all trajectories. @type
            % rowvec<double> @default [1e-1 1e-2 1e-3 1e-4]
            % orders: The different orders for which to compute the minimum
            % ErrorOrders @type rowvec<integer> @default One to
            % approx.DEIM.MaxOrder in steps of three
            % mu: The parameters that determine the desired trajectories.
            % The minimum needed M' will be averaged over all trajectories.
            % @type matrix<double> @default All full model's parameter
            % samples
            %
            % Return values:
            % res: A matrix containing the minimum required M', indexed by
            % orders in rows and relative errors in columns. @type
            % matrix<double>
            % relerrs: The effectively used relerrs. @type rowvec<double>
            % orders: The effectively used orders. @type rowvec<integer>
            % t: A PrintTable containing the results. If not specified as a
            % nargout argument, the table will be printed instead. @type
            % PrintTable @optional
            fm = r.FullModel;
            f = fm.Approx;
            if nargin < 5
                inidx = [];
            end
            if nargin < 4 || isempty(mu)
                mu = fm.Data.ParamSamples(:,1:10);
            end
            if nargin < 3 || isempty(orders)
                orders = 1:3:f.MaxOrder;
            end
            if nargin < 2 || isempty(relerrs)
                relerrs = [1e-1 1e-2 1e-3 1e-4];
            end
            
            oldo = f.Order;
            res = zeros(length(orders),length(relerrs));
            t = PrintTable('Mean required M'', averaged over %d trajectories',size(mu,2));
            t.HasHeader = true;
            title = arrayfun(@(e)sprintf('%g',e),relerrs,'Unif',false);
            t.addRow('DEIM order / rel. error',title{:});
            pi = tools.ProcessIndicator('Computing mean required error orders for %d DEIM orders and %d relative errors',...
                length(orders)*length(relerrs),...
                false,length(orders),length(relerrs));
            for i = 1:length(orders)
                f.Order = orders(i);
                [etrue, ~, ED] = testing.DEIM.getTrajApproxErrorDEIMEstimates(r, mu, inidx);
                ET = repmat(etrue',[1 1 size(ED,1)]);
                ET = shiftdim(ET,2);
                if size(ET,2) == 1 && size(ET,1) > 1
                    ET = ET';
                end
                if size(ED,1) == 1
                    ED = squeeze(ED);
                end
                rel = max(ED ./ ET,[],2);
                for j = 1:length(relerrs)
                    % Trick: find nonzero flags, get first occurrence for each
                    % column and use their row index as index of the ErrorOrder
                    % at which the relative error (between estimate and true)
                    % is smaller than relerrs(i) for all time-steps!
                    [v, c] = find(rel < relerrs(j));
                    [~, firstpos] = unique(c,'first');
                    res(i,j) = mean(v(firstpos));
                    pi.step;
                end
                tmp = res(i,:);
                tmp(isnan(tmp)) = -1;
                res(i,:) = tmp;
                hlp = num2cell(res(i,:));
                t.addRow(orders(i),hlp{:},{'%g'});
            end
            pi.stop;
            f.Order = oldo;
            
            if nargout < 4
                t.display;
            end
        end
        
%         function [res, relerrs, orders, t] = getMeanRequiredErrorOrders(r, relerrs, orders)
%             % [res, relerrs, orders, t] = getMeanRequiredErrorOrders(this, relerrs, orders, mu)
%             % This high-level function computes the minimum values of M' in
%             % order to have the relative error of true to estimated DEIM
%             % error smaller than the specified relerrs.
%             %
%             % Parameters:
%             % relerrs: The relative errors that may occur between true and
%             % estimated approximation error, maximized over time of each
%             % trajectory, averaged over all trajectories. @type
%             % rowvec<double> @default [1e-1 1e-2 1e-3 1e-4]
%             % orders: The different orders for which to compute the minimum
%             % ErrorOrders @type rowvec<integer> @default One to
%             % approx.DEIM.MaxOrder in steps of three
%             % mu: The parameters that determine the desired trajectories.
%             % The minimum needed M' will be averaged over all trajectories.
%             % @type matrix<double> @default All full model's parameter
%             % samples
%             %
%             % Return values:
%             % res: A matrix containing the minimum required M', indexed by
%             % orders in rows and relative errors in columns. @type
%             % matrix<double>
%             % relerrs: The effectively used relerrs. @type rowvec<double>
%             % orders: The effectively used orders. @type rowvec<integer>
%             % t: A PrintTable containing the results. If not specified as a
%             % nargout argument, the table will be printed instead. @type
%             % PrintTable @optional
%             fm = r.FullModel;
%             f = fm.Approx;
%             if nargin < 3 || isempty(orders)
%                 orders = 1:3:f.MaxOrder;
%             end
%             if nargin < 2 || isempty(relerrs)
%                 relerrs = [1e-1 1e-2 1e-3 1e-4];
%             end
%             
%             oldo = f.Order;
%             res = zeros(length(orders),length(relerrs));
%             t = PrintTable('Mean required M'', averaged over %d trajectories',size(mu,2));
%             t.HasHeader = true;
%             title = arrayfun(@(e)sprintf('%g',e),relerrs,'Unif',false);
%             t.addRow('DEIM order / rel. error',title{:});
%             pi = tools.ProcessIndicator('Computing mean required error orders for %d DEIM orders and %d relative errors',...
%                 length(orders)*length(relerrs),...
%                 false,length(orders),length(relerrs));
%             for i = 1:length(orders)
%                 f.Order = orders(i);
%                 [etrue, ~, ED] = testing.DEIM.getTrajApproxErrorDEIMEstimates(r, mu, inidx);
%                 ET = repmat(etrue',[1 1 size(ED,1)]);
%                 ET = shiftdim(ET,2);
%                 if size(ET,2) == 1 && size(ET,1) > 1
%                     ET = ET';
%                 end
%                 if size(ED,1) == 1
%                     ED = squeeze(ED);
%                 end
%                 rel = max(ED ./ ET,[],2);
%                 for j = 1:length(relerrs)
%                     % Trick: find nonzero flags, get first occurrence for each
%                     % column and use their row index as index of the ErrorOrder
%                     % at which the relative error (between estimate and true)
%                     % is smaller than relerrs(i) for all time-steps!
%                     [v, c] = find(rel < relerrs(j));
%                     [~, firstpos] = unique(c,'first');
%                     res(i,j) = mean(v(firstpos));
%                     pi.step;
%                 end
%                 tmp = res(i,:);
%                 tmp(isnan(tmp)) = -1;
%                 res(i,:) = tmp;
%                 hlp = num2cell(res(i,:));
%                 t.addRow(orders(i),hlp{:},{'%g'});
%             end
%             pi.stop;
%             f.Order = oldo;
%             
%             if nargout < 4
%                 t.display;
%             end
%         end
        
        %% Model DEIM reduction quality assessment pics
        function [errs, relerrs, times, deim_orders] = getDEIMReducedModelErrors(r, mu, inidx, deim_orders)
            d = r.System.f;
            if nargin < 4
                deim_orders = 1:d.MaxOrder;
            end
            oldo = d.Order;
            olde = r.ErrorEstimator.Enabled;
            no = length(deim_orders);
            errs = zeros(no,length(r.Times));
            relerrs = errs;
            times = zeros(no+1,1);
            [~, y, ct] = r.FullModel.simulate(mu, inidx);
            times(end) = ct;
            r.ErrorEstimator.Enabled = false;
            pi = tools.ProcessIndicator('Computing DEIM-reduced model simulation errors for %d orders',no,false,no);
            for m = 1:no
                r.System.f.Order = deim_orders(m);
                [~, yr, ct] = r.simulate(mu, inidx);
                errs(m,:) = Norm.L2(y-yr);
                relerrs(m,:) = errs(m,:)./Norm.L2(y);
                times(m) = ct;
                pi.step;
            end
            pi.stop;
            d.Order = oldo;
            r.ErrorEstimator.Enabled = olde;
        end
        
        function pm = getDEIMReducedModelErrors_plots(r, errs, relerrs, times, deim_orders, pm)
            if nargin < 6
                pm = tools.PlotManager(false);
            end
            [X, Y] = meshgrid(r.Times, deim_orders);
            h = pm.nextPlot('abserr',sprintf(['L2-absolute reduction errors\n'...
                '(Linf in time for original view)']),...
                'time','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,errs,'EdgeColor','interp');
            view(90,0);
            
            h = pm.nextPlot('relerr','L2-relative reduction errors',...
                'time','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,relerrs,'EdgeColor','interp');
            view(-120,30);
            
            h = pm.nextPlot('ctimes','Computation times for reduced models',...
                'DEIM order');
            plot(h,deim_orders,times(1:end-1));
            
            h = pm.nextPlot('speedup',...
                sprintf('Speedup against full model simulation of %gs',times(end)),...
                'DEIM order');
            plot(h,deim_orders,times(end)./times(1:end-1));
            
            if nargout == 0
                pm.done;
            end
        end
        
        %% Matrix DEIM approximation analysis
        function [e, aln, orders] = compareDEIM_Full_Jacobian(m, atd, orders)
            deim = m.ErrorEstimator.JacMDEIM;
            if nargin < 3
                orders = 1:deim.MaxOrder;
                if nargin < 2
                    atd = m.Data.JacobianTrainData;
                end
            end
            mo = length(orders);
            n = size(atd.xi,2);
            e = zeros(mo,n,3);
            aln = zeros(mo,n);
            oldo = deim.Order;
            Qk = deim.Qk;
            deim.setSimilarityTransform([]);
            pi = tools.ProcessIndicator('Computing matrix DEIM approximation error and log norm errors for %d orders',mo,false,mo);
            JP = logical(m.System.f.JSparsityPattern);
            for o = 1:mo
                deim.Order = o;
                for i=1:n
                    ti = atd.ti(i);  mui = atd.mui(:,i); 
                    J = m.System.f.getStateJacobian(atd.xi(:,i),ti,mui);
                    DJ = reshape(deim.evaluate(m.Data.V'*atd.xi(:,i),ti,mui),...
                        size(J,1),[]); 
                    diff = J-DJ;
                    diff = diff(JP);
                    e(o,i,1) = max(max(abs(diff)));
                    e(o,i,2) = Norm.L2(diff);
                    e(o,i,3) = mean(abs(diff));
                    e(o,i,4) = var(diff);

                    % Get log norm
                    aln(o,i) = general.Utils.logNorm(DJ);
                end
                pi.step;
            end
            pi.stop;
            deim.Order = oldo;
            deim.setSimilarityTransform(Qk);
        end
        
        function pm = compareDEIM_Full_Jacobian_plots(m, e, aln, orders, pm, atdsubset)
            if nargin < 5 || isempty(pm)
                pm = tools.PlotManager(false,1,2);
            end
            co = 'none';
            td = 1:size(e,2);
            [X,Y] = meshgrid(td,orders);
            e = sort(e,2);
            h = pm.nextPlot('max_diff','Maximum absolute difference',...
                'trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,e(:,:,1),'EdgeColor',co);
            h = pm.nextPlot('max_diff_mean','Mean maximum absolute difference','DEIM order');
            semilogy(h,orders,mean(e(:,:,1),2));
            
            h = pm.nextPlot('l2_diff','L2 vector difference',...
                'trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,e(:,:,2),'EdgeColor',co);
            h = pm.nextPlot('l2_diff_mean','Mean L2 vector difference','DEIM order');
            semilogy(h,orders,mean(e(:,:,2),2));
            
            h = pm.nextPlot('mean_diff','Mean difference value',...
                'trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,e(:,:,3),'EdgeColor',co);
            h = pm.nextPlot('mean_diff_mean','Mean mean difference value','DEIM order');
            semilogy(h,orders,mean(e(:,:,3),2));
            
            h = pm.nextPlot('var_diff','Variance of difference',...
                'trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,e(:,:,4),'EdgeColor',co);
            h = pm.nextPlot('var_diff_mean','Mean variance of difference','DEIM order');
            semilogy(h,orders,mean(e(:,:,4),2));
            
            ln = m.Data.JacSimTransData.LogNorms;
            if length(ln) ~= size(aln,2)
                ln = ln(atdsubset);
            end
            ln = repmat(ln,length(orders),1);
            err = sort(abs(ln-aln),2);
            h = pm.nextPlot('log_norm_diff','Differences of logarithmic norms',...
                'trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,err,'EdgeColor','none');
            
            err = sort(abs((ln-aln) ./ ln),2);
            h = pm.nextPlot('log_norm_diff','Differences of logarithmic norms',...
                'trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,err,'EdgeColor','none');
            
            if nargout < 1
                pm.done;
            end
        end
        
        %% Comparison of similarity transformed jacobian log norms to full log norms
        function [aln, times, st_sizes] = compareSimTransJac_FullJac(m, st_sizes)
            % Computes logarithmic norms of similarity transformed
            % jacobians using the model's offline data, containing `N`
            % trajectory samples.
            %
            % Parameters:
            % m: A BaseFullModel whos offlineGenerations have been run and
            % that uses a DEIMEstimator. @type models.BaseFullModel
            % st_sizes: The sizes `s_1,\ldots,s_n` to use for the
            % similarity transformation. If left empty, all from one to the
            % models ErrorEstimator.JacSimTransMaxSize are used
            % (extensive!) @type rowvec<integer>
            %
            % Return values:
            % aln: A `n \times N` matrix with the approximated logarithmic
            % norms in rows for each sim. trans. size.
            % times: A `n \times N` matrix with the computation time for
            % the logarithmic norm
            % st_sizes: The effectively used sizes (If st_sizes is given as
            % parameter, it is looped through)
            
            % Get required data
            jtd = m.Data.JacobianTrainData;
            QFull = m.ErrorEstimator.QFull;
            f = m.System.f;
            
            % Input checks
            if nargin < 2
                st_sizes = 1:size(QFull,2);
            end
            
            % Preps
            stn = length(st_sizes);
            n = size(jtd.xi,2);
            aln = zeros(stn,n);
            times = zeros(stn,n);
            pi = tools.ProcessIndicator('Computing approximated log norms for %d full jacobians and %d sim. trans. orders',...
                n*stn,false,n,stn);
            for nr = 1:n
                J = f.getStateJacobian(jtd.xi(:,nr),jtd.ti(nr),jtd.mui(:,nr));
                for snr = 1:stn
                    Qk = QFull(:,1:st_sizes(snr));
                    t = tic;
                    aln(snr,nr) = general.Utils.logNorm(Qk'*J*Qk);
                    times(snr,nr) = toc(t);
                    pi.step;
                end
            end
            pi.stop;
            if nargout == 0
                testing.DEIM.compareSimTransJac_FullJac_plots(m, aln, times, st_sizes);
            end
        end
        
        function pm = compareSimTransJac_FullJac_plots(m, aln, times, st_sizes, pm)
            if nargin < 5
                pm = tools.PlotManager(false,2,2);
            end
            jstd = m.Data.JacSimTransData;
            pm.nextPlot('correct_ln',sprintf('Full logarithmic norms\nJacobian full dimension %dx%d',...
                m.System.f.fDim,m.System.f.xDim));
            plot(jstd.LogNorms);
            
            if size(aln,1) > 1
                [X,Y] = meshgrid(1:size(aln,2),st_sizes);
                LN = repmat(jstd.LogNorms,size(aln,1),1);
            end
            
            if size(aln,1) == 1
                h = pm.nextPlot('approx_ln',...
                    sprintf('Approx logarithmic norms\nSimilarity transformation size %d',...
                    st_sizes));
                plot(h, aln);
            else
                h = pm.nextPlot('approx_ln','Approx logarithmic norms',...
                    'training sample','similarity transform size');
                tools.LogPlot.nicesurf(h, X, Y, aln);
            end
            
            h = pm.nextPlot('rel_err','Relative error','training sample',...
                'similarity transform size');
            tools.LogPlot.logsurf(h,X,Y,abs((LN-aln)./LN));
            
            sv = m.ErrorEstimator.QSingVals;
            s = m.ErrorEstimator.JacSimTransMaxSize;
            h = pm.nextPlot('sing_vals',...
                sprintf('Singular values of SVD of eigenvector matrix\nBlack line: Max size of sim. trans., singular value= %g',...
                sv(s)),'POD Modes','Singular values');
            semilogy(h,sv);
            hold on;
            plot(h,[s s+eps],[min(sv) max(sv)],'k');
            hold off;
            
            h = pm.nextPlot('comp_times',...
                sprintf('Average computation times (over %d values)\nfor log norms of similarity transformed matrices',...
                    size(times,2)),'similarity transformation size','mean computation time');
            me = mean(times,2);
            semilogy(h,st_sizes,me,'-s');
            if nargout < 1
                pm.done;
            end
        end
        
        %% Comparison of similarity transformed DEIM-approximated jacobian log norms to full log norms
        function [aln, times, jtimes, deim_orders, st_sizes] = ...
                compareSimTransDEIMJac_FullJac(m, deim_orders, st_sizes)
            % Computes logarithmic norms of similarity transformed AND
            % matrix DEIM approximated jacobians using the model's offline
            % data, containing `N` trajectory samples.
            %
            % Parameters:
            % m: A BaseFullModel whos offlineGenerations have been run and
            % that uses a DEIMEstimator. @type models.BaseFullModel
            % deim_orders: The DEIM orders `d_1,\ldots,d_m` to set for the
            % matrix DEIM of the jacobian. @type rowvec<integer>
            % st_sizes: The sizes `s_1,\ldots,s_n` to use for the
            % similarity transformation. If left empty, all from one to the
            % models ErrorEstimator.JacSimTransMaxSize are used
            % (extensive!) @type rowvec<integer>
            %
            % Return values:
            % aln: A `m \times n \times N` matrix with the approximated logarithmic
            % norms in rows for each sim. trans. size.
            % times: A `m \times n \times N` matrix with the computation times for
            % the logarithmic norm
            % jtimes: A `m \times n \times N` matrix with the computation
            % times for the sim. trans. DEIM approximated jacobians
            % deim_orders: The effectively used DEIM orders (If given as
            % parameter, it is looped through)
            % st_sizes: The effectively used sizes (If given as
            % parameter, it is looped through)
            
            % Get required data
            jtd = m.Data.JacobianTrainData;
            zi = m.Data.W'*jtd.xi;
            e = m.ErrorEstimator;
            jd = e.JacMDEIM;
            
            % Input checks
            if nargin < 3
                st_sizes = [1:e.JacSimTransMaxSize 0];
                if nargin < 2
                    deim_orders = 1:jd.MaxOrder;
                end
            end
            
            % Preps
            stn = length(st_sizes);
            n = size(jtd.xi,2);
            no = length(deim_orders);
            
            aln = zeros(no,stn,n);
            times = aln;
            jtimes = aln;
            pi = tools.ProcessIndicator('Computing approximated log norms over %d DEIM orders and %d sim. trans. sizes on %d training values',...
                numel(aln),false,no,stn,n);
            for onr = 1:no
                jd.Order = deim_orders(onr);
                for snr = 1:stn
                    jd.setSimilarityTransform(e.QFull(:,1:st_sizes(snr)));
                    for nr = 1:n
                        t = tic;
                        J = jd.evaluate(zi(:,nr),jtd.ti(nr),jtd.mui(:,nr));
                        jtimes(onr,snr,nr) = toc(t);
                        t = tic;
                        aln(onr,snr,nr) = general.Utils.logNorm(J);
                        times(onr,snr,nr) = toc(t);
                        pi.step;
                    end
                end
            end
            pi.stop;
            if nargout == 0
                testing.DEIM.compareSimTransDEIMJac_FullJac_plots(m, aln, times, st_sizes);
            end
        end
        
        function pm = compareSimTransDEIMJac_FullJac_plots(m, aln, ...
                times, jtimes, deim_orders, st_sizes, pm)
            if nargin < 7
                pm = tools.PlotManager(false,2,2);
            end
            
            jstd = m.Data.JacSimTransData;
            n = size(aln,3);
            ln = reshape(jstd.LogNorms,1,1,[]);
            ln = repmat(ln,[length(deim_orders), length(st_sizes), 1]);
            
            abserr = abs(aln-ln);
            relerr = abs(abserr ./ ln);
            relerr(isnan(relerr)) = 0;
            
            % Aggregate errors over trajectory sample data
            %aln = sqrt(sum(aln.^2,3));
            abserr = mean(abserr,3);
            relerr = mean(relerr,3);
            
            % Mean the times too
            ttimes = times + jtimes;
            times = mean(times,3);
            jtimes = mean(jtimes,3);
            ttimes = mean(ttimes,3);
            
            [X,Y] = meshgrid(st_sizes,deim_orders);
            
            h = pm.nextPlot('abs_err',sprintf('Mean absolute approximation error over %d samples',n),...
                'Similarity transformation size','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,abserr);
            
            h = pm.nextPlot('rel_err',sprintf('Mean relative approximation error over %d samples',n),...
                'Similarity transformation size','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,relerr);
            
            h = pm.nextPlot('comp_times_j',...
                sprintf('Average computation times over %d values\nfor matrix DEIM jacobian evaluation',...
                    n),'Similarity transformation size','DEIM order');
            tools.LogPlot.nicesurf(h,X,Y,jtimes);
            
            h = pm.nextPlot('comp_times',...
                sprintf('Average computation times over %d values\nfor log norm computation of sim.trans. matrix DEIM jacobian',...
                    n),'Similarity transformation size','DEIM order');
            tools.LogPlot.nicesurf(h,X,Y,times);
            
            h = pm.nextPlot('comp_times_total',...
                sprintf('Average total computation times over %d values\n(jac comp + log norm comp)',...
                    n),'Similarity transformation size','DEIM order');
            tools.LogPlot.nicesurf(h,X,Y,ttimes);
            
            if nargout < 1
                pm.done;
            end
        end
        
        function [m, r] = test_DEIM
            m = models.pcd.PCDModel(1);
            m.EnableTrajectoryCaching = false;
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 40;
            m.System.Params(1).Desired = 10;
            m.SpaceReducer = spacereduction.PODGreedy;
            m.offlineGenerations;
            m.Approx.Order = [20 4];
            r = m.buildReducedModel;
            save approx.test_DEIM m r;
        end
    end    
end
