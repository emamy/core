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
            ma = ModelAnalyzer(m.buildReducedModel);
            
            mu = m.Data.ParamSamples(:,5);
%             [t, pm] = ma.compareRedFull(mu);
%             t.Format = 'tex';
%             t.saveToFile('comp_red_full.tex');
%             pm.savePlots('.',{'fig','jpg','eps'}, true);
            
            d = m.Approx;
            res = testing.DEIM.computeDEIMErrors(d, m.Data.ApproxTrainData);
            %pm = PlotManager(true);
            pm = PlotManager(false,1,2);
            pm.FilePrefix = 'DEIM_errors';
            testing.DEIM.plotDEIMErrs(res, pm);
            %pm.savePlots('.',{'fig','jpg','eps'}, true);
            pm.savePlots('.',{'fig','jpg'}, true);
            
            [etrue, EE, ED] = ma.getTrajApproxErrorDEIMEstimates(mu,[]);
%             pm = PlotManager(true);
            pm = PlotManager(false,1,2);
            pm.FilePrefix = 'traj_approx_err';
            ma.getTrajApproxErrorDEIMEstimates_plots(m.scaledTimes, etrue, EE, ED, pm);
%             pm.savePlots('.',{'fig','jpg','eps'}, true);
            pm.savePlots('.',{'fig','jpg'}, true);
            
            [minreq, relerrs, orders, t] = ma.getMeanRequiredErrorOrders;
            t.Format = 'tex';
            t.saveToFile('meanrequirederrororders.tex');
            
            save analysis_DEIM_approx;
        end
        
        %% DEIM approximation analysis over training data
        function res = computeDEIMErrors(deim, atd, orders, errorders)
            % Computes the DEIM approximation error over the given training
            % data for specified DEIM orders (M) and DEIM error orders
            % (M').
            %
            % The error over any training data point is measured L2 in
            % state and Linf over all training points.
            %
            % Parameters:
            % deim: The DEIM instance @type general.DEIM
            % atd: The approximation training data @type data.ApproxTrainData
            % orders: The different orders `1 \leq m\leq M` to try @type rowvec<integer>
            % errorders: The different error orders `1 \leq m'\leq M-m` to try @type
            % rowvec<integer>
            %
            % Return values:
            % res: A 12 x totalcombinations matrix containing the values in
            % rows
            % 1: order
            % 2: errororder
            % 3: max abs true error (indep. of errororder)
            % 4: max rel true error w.r.t true fxi norm (indep. of errororder)
            % 5: mean abs true error (indep. of errororder)
            % 6: mean rel true error w.r.t true fxi norm (indep. of errororder)
            % 7: max abs estimated error
            % 8: max rel estimated error w.r.t true fxi norm
            % 9: mean abs estimated error
            % 10: mean rel estimated error w.r.t true fxi norm
            % 11: max rel error of (estimated-true error) w.r.t true error
            % 12: mean rel error of (estimated-true error) w.r.t true error
            % 13: max rel error of (estimated-maxorder error) w.r.t maxorder error
            % 14: mean rel error of (estimated-maxorder error) w.r.t maxorder error
            oldo = deim.Order;
            if nargin < 4
                if nargin < 3
                    orders = 1:deim.MaxOrder-1;
                end
                errorders = [];
                neworders = [];
                for i=1:length(orders)
                    new = 1:(deim.MaxOrder-orders(i));
                    errorders = [errorders new];%#ok
                    neworders = [neworders orders(i)*ones(1,length(new))];%#ok
                end
                orders = neworders;
            elseif ~isempty(errorders) && length(orders) ~= length(errorders)
                error('If both the orders and error orders are given they must have same number of elements');
            end
            n = length(orders);
            res = zeros(14,n);
            res(1,:) = orders;
            res(2,:) = errorders;
            % State space error function
            efun = @Norm.L2;
            
            nb = atd.xi.nBlocks;
            pi = ProcessIndicator(sprintf('Computing DEIM errors and estimates for %d Order/ErrOrder settings over %d xi-Blocks',n,nb),n*nb);
            co = [];
            m = atd.xi.m;
            for k = 1:nb
                pos = atd.xi.getBlockPos(k);
                xi = atd.xi(:,pos);
                fxi = atd.fxi(:,pos); % must have same length, but may have different block size.
                if ~isempty(deim.V)
                    xi = deim.V'*xi;
                end
                fxinorm = efun(fxi);
                deim.Order = deim.MaxOrder;
                afximaxorder = deim.evaluate(xi,atd.ti(pos),atd.mui(:,pos));
                for i = 1:n
                    deim.Order = [res(1,i) res(2,i)];
                    % Check if new overall DEIM order m
                    if isempty(co) || res(1,i) ~= co
                        afxi = deim.evaluate(xi,atd.ti(pos),atd.mui(:,pos));
                        % Get true DEIM error against full function
                        etrue = efun(fxi-afxi);
                        % Get error measured against max order DEIM approx
                        emaxorder = efun(afximaxorder-afxi);
                        % 3: max abs true error (indep. of errororder)
                        res(3,i) = max(res(3,i),max(etrue)); 
                        % 4: max rel true error w.r.t true fxi norm (indep. of errororder)
                        hlp = etrue./fxinorm;
                        hlp(isnan(hlp)) = 0;
                        res(4,i) = max(res(4,i),max(hlp)); 
                        % 5: mean abs true error (indep. of errororder)
                        res(5,i) = res(5,i) + sum(etrue)/m; 
                        % 6: mean rel true error w.r.t true fxi norm (indep. of errororder)
                        res(6,i) = res(6,i) + sum(hlp)/m; 
                        co = res(1,i);
                    else % else copy
                        res(3:6,i) = res(3:6,i-1);
                    end
                    if res(2,i) > 0
                        % Estimated absolute/rel errors
                        eest = efun(deim.getEstimatedError(xi,atd.ti(pos),atd.mui(:,pos)));
                        % 7: max abs estimated error
                        res(7,i) = max(res(7,i),max(eest));
                        % 8: max rel estimated error w.r.t true fxi norm
                        hlp = eest./fxinorm;
                        hlp(isnan(hlp)) = 0;
                        res(8,i) = max(res(8,i),max(hlp));
                        % 9: mean abs estimated error
                        res(9,i) = res(9,i) + sum(eest)/m;
                        % 10: mean rel estimated error w.r.t true fxi norm
                        res(10,i) = res(10,i) + sum(hlp)/m;
                        % 11: max rel error on (estimated-true error) w.r.t true error
                        hlp = abs(eest-etrue)./etrue;
                        hlp(isnan(hlp)) = 0;
                        res(11,i) = max(res(11,i),max(hlp));
                        % 12: mean rel error on (estimated-true error) w.r.t true error
                        res(12,i) = res(12,i) + sum(hlp)/m;
                        % 13: max rel error on (estimated-maxorder error) w.r.t maxorder error
                        hlp = abs(eest-emaxorder)./emaxorder;
                        hlp(isnan(hlp)) = 0;
                        res(13,i) = max(res(13,i),max(hlp));
                        % 14: mean rel error on (estimated-maxorder error) w.r.t maxorder error
                        res(14,i) = res(14,i) + sum(hlp)/m;
                    end
                    pi.step;
                end
            end
            deim.Order = oldo;
            pi.stop;
        end
        
        function pm = plotDEIMErrs(res, pm)
            
            if nargin < 2
                pm = PlotManager(false,2,3);
                pm.FilePrefix = 'comp_m_mdash';
            end
            % Not used yet:
            % 5: mean abs true error (indep. of errororder)
            % 6: mean rel true error w.r.t true fxi norm (indep. of errororder)
            % 9: mean abs estimated error
            % 10: mean rel estimated error w.r.t true fxi norm
            
            % Compute unique positions
            [orders, idx] = unique(res(1,:));
            % 1: order, 2: errororder
            tri = delaunay(res(1,:),res(2,:));
            
            % 3: max abs true error (indep. of errororder)
            h = pm.nextPlot('true_abs','Linf-L2 absolute error','m','m''');
            semilogy(h,orders,res(3,idx));
            %doplot(h,tri,res(1,:),res(2,:),res(3,:));
            % 7: max abs estimated error
            h = pm.nextPlot('est_abs','Estimated absolute error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),res(7,:),'EdgeColor','none');
            
            % 4: max rel true error w.r.t true fxi norm (indep. of errororder)
            h = pm.nextPlot('true_rel','Linf-L2 relative  error','m','m''');
            semilogy(h,orders,res(4,idx));
            %doplot(h,tri,res(1,:),res(2,:),res(4,:));
            % 8: max rel estimated error w.r.t true fxi norm
            h = pm.nextPlot('est_rel','Estimated relative error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),res(8,:),'EdgeColor','none');
            
            %% Reduced Triplots
            n = size(res,2);
            %sel = 1:step:n;
            %sel = round(logspace(log10(1),log10(n),300));
            sel = [1 2 3 4 5 7 9 11 13 15:3:30 31:20:n];
            res = res(:,sel);
            n = size(res,2);
            % 1: order
            % 2: errororder
            tri = delaunay(res(1,:),res(2,:));
            
            % 11: max rel error of (estimated-true error) w.r.t true error
            % 12: mean rel error of (estimated-true error) w.r.t true error
            h = pm.nextPlot('max_relerr_est_true_totrue','"max (true - estimated)/true" relative error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),res(11,:));
            hold on;
            sh = doplot(h,tri,res(1,:),res(2,:),1e-2*ones(1,n),...
                'EdgeColor','none','FaceColor','k');
            alpha(sh,.2);
            hold off;
            axis ij;
            h = pm.nextPlot('mean_relerr_est_true_totrue','"mean (true - estimated)/true" relative error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),res(12,:));
            hold on;
            sh = doplot(h,tri,res(1,:),res(2,:),1e-2*ones(1,n),...
                'EdgeColor','none','FaceColor','k');
            alpha(sh,.2);
            hold off;
            axis ij;
            view(-40,34);
            
            % 13: max rel error of (estimated-maxorder error) w.r.t maxorder error
            % 14: mean rel error of (estimated-maxorder error) w.r.t maxorder error
            h = pm.nextPlot('max_relerr_est_emaxo_toemaxo','"max (maxordererror - estimated)/maxordererror" relative error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),res(13,:));
            hold on;
            sh = doplot(h,tri,res(1,:),res(2,:),1e-2*ones(1,n),...
                'EdgeColor','none','FaceColor','k');
            alpha(sh,.2);
            hold off;
            axis ij;
            h = pm.nextPlot('mean_relerr_est_emaxo_toemaxo','"mean (maxordererror - estimated)/maxordererror" relative error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),res(14,:));
            hold on;
            sh = doplot(h,tri,res(1,:),res(2,:),1e-2*ones(1,n),...
                'EdgeColor','none','FaceColor','k');
            alpha(sh,.2);
            hold off;
            axis ij;
            view(-40,34);
            
            h = pm.nextPlot('abs_diff_est_true','Absolute difference of estimated to true error','m','m''');
            doplot(h,tri,res(1,:),res(2,:),abs(res(3,:)-res(7,:)));
            view(75,34);

            pm.done;
            
            function p = doplot(h, tri, x, y, z, varargin)
                p = LogPlot.logtrisurf(h, tri, x, y, z, varargin{:});
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
                pm = PlotManager(false,3,1);
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
        
        function [efull, ered, fxno] = getApproxErrorFullRed(r, xr, t, mu, V)
            fm = r.FullModel;
            fm.Approx.Order = r.System.f.Order;
            fx = fm.System.f.evaluate(V*xr,t,mu);
            efull = Norm.L2(fx-r.FullModel.Approx.evaluate(V*xr,t,mu));
            ered = Norm.L2(fx-V*r.System.f.evaluate(xr,t,mu));
            fxno = Norm.L2(fx);
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
            % Parameters:
            % r: The reduced model @type models.ReducedModel
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
            % requested as output, no plotting will be done.
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
            
            ma = ModelAnalyzer(r);
            
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
                pm = testing.DEIM.getTrajApproxErrorDEIMEstimates_plots(r, etrue(1,:), EE(:,:,1), ED(:,:,1));
            end
        end
        
        function pm = getTrajApproxErrorDEIMEstimates_plots(r, etrue, EE, ED, pm)
            % Visualizes the results of the
            % ModelAnalyzer.getTrajApproxErrorDEIMEstimates method.
            f = r.FullModel.Approx;
            o = f.Order(1);
            mo = f.MaxOrder;
            if nargin < 6
                pm = PlotManager(false,2,2);
            end
            t = r.Times;
            h = pm.nextPlot('true_err','True approximation error on trajectory','time','error');
            semilogy(h,t,etrue);
            h = pm.nextPlot('est_err',sprintf('Estimated approximation errors on trajectory\nDEIM order = %d / max %d',f.Order(1),f.MaxOrder),...
                'time', 'DEIM error order');
            [T,O] = meshgrid(t,1:mo-o);
            LogPlot.logsurf(h,T,O,EE,'EdgeColor','none');
            view(-45,45);
            
            h = pm.nextPlot('abs_diff',sprintf('Difference true to estimated error\nDEIM order = %d / max %d',f.Order(1),f.MaxOrder),...
                'time','DEIM error order');
            LogPlot.logsurf(h,T,O,ED,'EdgeColor','none');
            view(-135,45);
            
            h = pm.nextPlot('rel_diff',sprintf('Relative error between true to estimated error\nDEIM order = %d / max %d',f.Order(1),f.MaxOrder),...
                'time','DEIM error order');
            LogPlot.logsurf(h,T,O,ED./repmat(etrue,mo-o,1),'EdgeColor','none','FaceColor','interp');
            hold on;
            p = surf(h,T,O,-ones(size(T))*2,'EdgeColor','red','FaceColor','k');
            alpha(p,.2);
            view(-135,45);
            hold off;
            pm.done;
        end
        
        function [t, values] = getMinReqErrorOrdersTable(errordata, relerrs, tsize, maxorder)
            % Return values:
            % t: A PrintTable containing the results. If not specified as a
            % nargout argument, the table will be printed instead. @type
            % PrintTable
            % values: a value matrix with rows for orders and columns for relerrs @type
            % matrix<double>
            if nargin == 4
                txt = sprintf('DEIM MaxOrder %d error',maxorder);
                maxvidx = 13;
                meanvidx = 14;
            else
                txt = 'true DEIM error';
                maxvidx = 11;
                meanvidx = 12;
            end
            t = PrintTable('Required M'' over %d training samples for min/mean errors relative to %s',...
                tsize,txt);
            t.HasHeader = true;
            title = arrayfun(@(e)sprintf('%g',e),relerrs,'Unif',false);
            t.addRow('m / rel. error',title{:});
            nr = length(relerrs);
            orders = unique(errordata(1,:));
            values = zeros(length(orders),1+2*nr);
            for i = 1:length(orders)
                o = orders(i);
                values(i,1) = o;
                pos = errordata(1,:) == o;
                maxv = errordata(maxvidx,pos);
                meanv = errordata(meanvidx,pos);
                hlp = cell.empty(0,nr);
                for j = 1:nr
                    % Max rel errs
                    maxidx = find(maxv <= relerrs(j),1,'first');
                    if isempty(maxidx), 
                        maxidx = min(maxv); 
                    end
                    % Mean rel errs
                    meanidx = find(meanv <= relerrs(j),1,'first');
                    if isempty(meanidx), 
                        meanidx = min(meanv); 
                    end
                    hlp{j} = sprintf('%g/%g',maxidx,meanidx);
                    values(i,1+[j j+nr]) = [maxidx meanidx];
                end
                t.addRow(o,hlp{:});
            end
            if nargout < 1
                t.display;
            end
        end
        
        %% DEIM approximation analysis over parameters for specific state space location
        function [mui, fxi, afxi] = getDEIMErrorsAtXForParams(m, x, numExtraSamples)
            %
            % Parameters:
            % m: The full model @type models.BaseFullModel
            % x: The state space location `\vx` @type colvec<double>
            % numExtraSamples: The number of extra samples to use. @type integer @default 0
            %
            % Only t=0 is used
            if nargin < 3
                numExtraSamples = 0;
            end
            mui = m.Data.ParamSamples;
            s = sampling.RandomSampler;
            s.Seed = 1;
            s.Samples = numExtraSamples;
            mui = [mui s.generateSamples(m)];
            xi = repmat(x,1,size(mui,2));
            ti = zeros(1,size(mui,2));
            fxi = m.System.f.evaluate(xi,ti,mui);
            afxi = m.Approx.evaluate(xi,ti,mui);
        end
        
        function pm = getDEIMErrorsAtXForParams_plots(m, mui, fxi, afxi, pm)
            if nargin < 5
                pm = PlotManager;
                pm.LeaveOpen = true;
            end

            tri = delaunay(mui(1,:),mui(2,:));
            
            err = Norm.L2(fxi-afxi);
            h = pm.nextPlot('abserr','Absolute errors over mu range');
            %trisurf(tri,mui(1,:),mui(2,:),err,'Parent',h);
            LogPlot.logtrisurf(h,tri,mui(1,:),mui(2,:),err);
            n = m.Data.SampleCount;
            hold on;
            plot3(mui(1,1:n),mui(2,1:n),log10(err(1:n)),'rx','MarkerSize',16);
            view(-8,-20);
            hold off;
            
            err = Norm.L2(fxi-afxi)./Norm.L2(fxi);
            h = pm.nextPlot('relerr','Relative errors over mu range');
            %trisurf(tri,mui(1,:),mui(2,:),err,'Parent',h);
            LogPlot.logtrisurf(h,tri,mui(1,:),mui(2,:),err);
            n = m.Data.SampleCount;
            hold on;
            plot3(mui(1,1:n),mui(2,1:n),log10(err(1:n)),'rx','MarkerSize',16);
            hold off;
            view(0,0);
        end
        
        %% Model DEIM reduction quality assessment pics        
        function [errs, relerrs, times, deim_orders] = getDEIMReducedModelErrors(r, mu, inidx, deim_orders)
            %
            % Parameters:
            % r: The reduced model @type models.ReducedModel
            % mu: The current parameter `\vmu` @type colvec<double>
            % inidx: The input index `i` for the input function `\vu_i(t)` @type integer
            % deim_orders: The DEIM orders `1 \leq m \leq M` to use @type rowvec<integer>
            d = r.System.f;
            if nargin < 4
                deim_orders = 1:(d.MaxOrder-r.System.f.Order(2));
            end
            oldo = d.Order;
            no = length(deim_orders);
            errs = zeros(no,length(r.Times));
            relerrs = errs;
            times = zeros(no+1,1);
            [~, y, ct] = r.FullModel.simulate(mu, inidx);
            times(end) = ct;
            pi = ProcessIndicator('Computing DEIM-reduced model simulation errors for %d orders',no,false,no);
            for m = 1:no
                r.System.f.Order = [deim_orders(m) r.System.f.Order(2)];
                [~, ~, ct] = r.simulate(mu, inidx);
                errs(m,:) = r.ErrorEstimator.OutputError;
                relerrs(m,:) = errs(m,:)./Norm.L2(y);
                times(m) = ct;
                pi.step;
            end
            pi.stop;
            d.Order = oldo;
        end
        
        function pm = getDEIMReducedModelErrors_plots(r, errs, relerrs, times, deim_orders, pm, tag)
            if nargin < 7
                tag = '';
                if nargin < 6
                    pm = PlotManager(false);
                end
            else
                ftag = [tag '_'];
            end
            [X, Y] = meshgrid(r.Times, deim_orders);
            h = pm.nextPlot([ftag 'abserr'],sprintf(['L2-absolute reduction errors\n'...
                '(Linf in time for original view), tag:' tag]),...
                'time','DEIM order');
            LogPlot.logsurf(h,X,Y,errs,'EdgeColor','interp');
            view(90,0);
            
            h = pm.nextPlot([ftag 'relerr'],['L2-relative reduction errors, tag:' tag],...
                'time','DEIM order');
            LogPlot.logsurf(h,X,Y,relerrs,'EdgeColor','interp');
            view(-120,30);
            
%             h = pm.nextPlot([ftag 'ctimes'],'Computation times for reduced models',...
%                 'DEIM order');
%             plot(h,deim_orders,times(1:end-1));
%             
%             h = pm.nextPlot([ftag 'speedup'],...
%                 sprintf('Speedup against full model simulation of %gs',times(end)),...
%                 'DEIM order');
%             plot(h,deim_orders,times(end)./times(1:end-1));
%             
%             if nargout == 0
%                 pm.done;
%             end
        end
        
        %% Matrix DEIM approximation analysis
        function [e, aln, orders] = compareDEIM_Full_Jacobian(m, atd, orders)
            % Compares the MDEIM approximation with the full jacobian
            %
            % Parameters:
            % m: The full model @type models.BaseFullModel
            % atd: The approximation training data @type data.ApproxTrainData
            % orders: The different orders `m_J` to try @type rowvec<integer>
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
            pi = ProcessIndicator('Computing matrix DEIM approximation error and log norm errors for %d orders',mo,false,mo);
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
                    aln(o,i) = Utils.logNorm(DJ);
                end
                pi.step;
            end
            pi.stop;
            deim.Order = oldo;
            deim.setSimilarityTransform(Qk);
        end
        
        function pm = compareDEIM_Full_Jacobian_plots(m, e, aln, orders, pm, atdsubset)
            if nargin < 5 || isempty(pm)
                pm = PlotManager(false,1,2);
            end
            co = 'none';
            td = 1:size(e,2);
            [X,Y] = meshgrid(td,orders);
            e = sort(e,2);
            h = pm.nextPlot('max_diff','Maximum absolute difference',...
                'trajectory point','DEIM order');
            LogPlot.logsurf(h,X,Y,e(:,:,1),'EdgeColor',co);
            h = pm.nextPlot('max_diff_mean','Mean maximum absolute difference','DEIM order');
            semilogy(h,orders,mean(e(:,:,1),2));
            
            h = pm.nextPlot('l2_diff','L2 vector difference',...
                'trajectory point','DEIM order');
            LogPlot.logsurf(h,X,Y,e(:,:,2),'EdgeColor',co);
            h = pm.nextPlot('l2_diff_mean','Mean L2 vector difference','DEIM order');
            semilogy(h,orders,mean(e(:,:,2),2));
            
            h = pm.nextPlot('mean_diff','Mean difference value',...
                'trajectory point','DEIM order');
            LogPlot.logsurf(h,X,Y,e(:,:,3),'EdgeColor',co);
            h = pm.nextPlot('mean_diff_mean','Mean mean difference value','DEIM order');
            semilogy(h,orders,mean(e(:,:,3),2));
            
            h = pm.nextPlot('var_diff','Variance of difference',...
                'trajectory point','DEIM order');
            LogPlot.logsurf(h,X,Y,e(:,:,4),'EdgeColor',co);
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
            LogPlot.logsurf(h,X,Y,err,'EdgeColor','none');
            
            err = sort(abs((ln-aln) ./ ln),2);
            h = pm.nextPlot('log_norm_diff','Differences of logarithmic norms',...
                'trajectory point','DEIM order');
            LogPlot.logsurf(h,X,Y,err,'EdgeColor','none');
            
            if nargout < 1
                pm.done;
            end
        end
        
        %% Effectivity analysis of error estimators
        function pm = effectivityAnalysis(r, mu, inputidx)
            % Plots an effectivity graph for the current error estimator.
            %
            % Parameters:
            % r: The reduced model @type models.ReducedModel
            pm = PlotManager(false,2,1);
            pm.SingleSize = [720 540];
            pm.LeaveOpen = true;
            [~,y] = r.FullModel.simulate(mu,inputidx);
            [t,yr] = r.simulate(mu,inputidx);
            err = Norm.L2(y-yr);
            h = pm.nextPlot('errors','True and estimated error','time','error');
            semilogy(h,t,err,'b',t,r.ErrorEstimator.OutputError,'r');
            h = pm.nextPlot('errors','True and estimated error','time','error');
            semilogy(h,t,r.ErrorEstimator.OutputError./err,'g');
        end
        
        %% Error estimator struct compilation
        function est = getDEIMEstimators_MDEIM_ST(rmodel, est, jdorders, stsizes)
            % Computes a set of DEIMEstimators for given MDEIM orders and simtrans-sizes.
            %
            % Parameters:
            % rmodel: The reduced model @type models.ReducedModel
            % est: The estimator settings struct @type struct
            % jdorders: The MDEIM orders `1\leq m_J \leq M_J` to use @type rowvec<integer>
            % stsizes: The similarity transformation sizes `1\leq k \leq d` to use @type
            % rowvec<integer>
            % 
            % Return values:
            % est: An estimator struct usable by the EstimatorAnalyzer @type struct
            if nargin < 3
                stsizes = [1 10];
                if nargin < 2
                    jdorders = [1 10];
                end
            end
            li = LineSpecIterator(length(jdorders)*length(stsizes)+4,1);
            for i = 1:length(est)
                li.excludeColor(est(i).Color);
            end
            
            % Different configurations
            e = rmodel.ErrorEstimator.clone;
            for j = 1:length(jdorders)
                cl = li.nextLineStyle;
                li2 = LineSpecIterator;
                for k = 1:length(stsizes)
                    str = sprintf('e.JacMDEIM.Order = %d; e.JacSimTransSize = %d;',...
                        jdorders(j),stsizes(k));
                    est(end+1).Name = sprintf('$m_J:%d, k:%d$',...
                        jdorders(j),stsizes(k));%#ok
                    est(end).Estimator = e;
                    est(end).Callback = @(e)eval(str);
                    est(end).MarkerStyle = li2.nextMarkerStyle;
                    est(end).LineStyle = cl;
                    est(end).Color = li.nextColor;
                end
            end
        end
        
        function est = getDEIMEstimators_ErrOrders(rmodel, est, errororders)
            if nargin < 2
                errororders = [1 2 5];
            end
            m = LineSpecIterator(2+length(errororders));
            for i = 1:length(est)
                m.excludeColor(est(i).Color);
            end
            % Different configurations
            e = rmodel.ErrorEstimator.clone; % need only one copy!
            for j = 1:length(errororders)
                str = sprintf('e.ReducedModel.System.f.Order = [e.ReducedModel.System.f.Order(1) %d];',errororders(j));
                est(end+1).Name = sprintf('$m''=%d$',errororders(j));%#ok
                est(end).Estimator = e;
                est(end).Callback = @(e)eval(str);
                est(end).MarkerStyle = m.nextMarkerStyle;
                est(end).LineStyle = '-.';
                est(end).Color = m.nextColor;
            end
        end
        
        %% General test methods
        function [res, m, r] = test_DEIM
            m = models.pcd.PCDModel(1);
            m.EnableTrajectoryCaching = false;
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 40;
            m.System.Params(1).Desired = 10;
            m.SpaceReducer = spacereduction.PODGreedy;
            m.offlineGenerations;
            m.Approx.Order = [20 4];
            r = m.buildReducedModel;
            %save approx.test_DEIM m r;
            res = true;
        end
    end    
end
