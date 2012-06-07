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
            res = d.computeDEIMErrors(m.Data.ApproxTrainData);
            %pm = tools.PlotManager(true);
            pm = tools.PlotManager(false,1,2);
            pm.FilePrefix = 'DEIM_errors';
            d.plotDEIMErrs(res, pm);
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
            oldo = deim.fOrder;
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
            
            res = zeros(6,no);
%             res = zeros(8,no);
            pi = tools.ProcessIndicator(sprintf('Computing DEIM errors and estimates for %d Order/ErrOrder settings',no),no);
            co = [];
            for i = 1:no
                o = orders(i);
                eo = errorders(i);
                deim.Order = [o eo];
                res(1:2,i) = [o; eo];
                
                % Absolute true error (only for new orders)
                if isempty(co) || o ~= co
                    afxi = deim.evaluate(atd.xi,atd.ti,atd.mui);
                    hlp = efun(atd.fxi-afxi);
                    res(3,i) = sumfun(hlp);
                    % Relative true error 
                    res(4,i) = sumfun(hlp./fxinorm);
                    co = o;
                else
                    res(3:4,i) = res(3:4,i-1);
                end
                % Estimated absolute/rel errors
                hlp = efun(deim.getEstimatedError(atd.xi,atd.ti,atd.mui));
                res(5,i) = sumfun(hlp);
                res(6,i) = sumfun(hlp./fxinorm);
                
%                 % Compute actual error between M and M' approximations
%                 deim.Order = [sum(deim.fOrder) 0];
%                 % Get order+errorder eval
%                 hlp = efun(deim.evaluateCoreFun(atd.xi,atd.ti,atd.mui) - afxi);
%                 res(7,i) = sumfun(hlp);
%                 res(8,i) = sumfun(hlp./fxinorm);
                pi.step;
            end
            pi.stop;
            if nargout == 2
                pm = testing.DEIM.plotDEIMErrs(res);
            end
            deim.Order = oldo;
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
               
        function [t, nof, nor, o, pm] = jacobian_tests(m, pm)
            % Tis function computes the jacobians of both the full and DEIM
            % approximated system functions.
            %
            % The maximum relative error for each DEIM order and
            % jacobian at a certain snapshot vector are computed and
            % displayed.
            atd = m.Data.ApproxTrainData;
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
        
        
        %% %%%%%%%%%%%%%%%%%%%% MATRIX DEIM STUFF %%%%%%%%%%%%%%%%%%%%
        function [no,no1,nom,jln,djln,sdjln] = matrix_deim(m, deim, nr)
            [x, mu] = m.Data.getTrajectoryNr(nr);
            t = m.scaledTimes;
            mu = repmat(mu,1,length(t));
            mo = deim.MaxOrder;
            VJ = m.Approx.VJ;
            pi = tools.ProcessIndicator('Computing stuff',mo*length(t));
            for o = 1:mo
                deim.Order = o;
                for i=1:length(t)
                    xi = x(:,i); 
                    ti = t(i); 
                    mui = mu(:,i); 
                    J = m.System.f.getStateJacobian(xi,ti,mui);
                    DJ = reshape(deim.evaluate(xi,ti,mui),size(J,1),[]); 
                    diff = J-DJ;
                    nz = J ~= 0;
                    no(o,i) = norm(full(diff)); %#ok<*AGROW>
                    no1(o,i) = norm(full(diff),1); %#ok<*AGROW>
                    rel = diff ./ J;
                    hlp = max(max(abs(rel(nz))));
                    if isempty(hlp)
                        hlp = 0;
                    end
                    nom(o,i) = hlp;
                    
                    jln(o,i) = general.Utils.logNorm(J);
                    djln(o,i) = general.Utils.logNorm(DJ);
                    sdjln(o,i) = general.Utils.logNorm(VJ'*DJ*VJ);
                    pi.step;
                end
            end
            pi.stop;
        end
        
        function pm = matrix_deim_plots(no,no1,nom,jln,djln,sdjln,pm)
            if nargin < 7
                pm = tools.PlotManager(false,1,3);
            end
            [X,Y] = meshgrid(1:size(no,2),1:size(no,1));
            h = pm.nextPlot('diff_norm','L2-Norm of difference matrix','trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,no);
            h = pm.nextPlot('diff_norm_l1','L1-Norm of difference matrix','trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,no1);
            h = pm.nextPlot('max_rel_diff','Maximum relative difference','trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,nom);
            
            h = pm.nextPlot('jac_log_norm',...
                'Log norm of full jacobian','trajectory point','DEIM order');
            surf(h,X,Y,jln,'FaceColor','interp','EdgeColor','interp');
            zlabel('log norm');
            h = pm.nextPlot('deim_jac_norm',...
                'Log norm of DEIM approximated jacobian',...
                'trajectory point','DEIM order');
            surf(h,X,Y,djln,'FaceColor','interp','EdgeColor','interp');
            zlabel('log norm');
            h = pm.nextPlot('st_deim_jac_norm',...
                sprintf('Log norm of partial similarity transformed\nDEIM approximated jacobian'),...
                'trajectory point','DEIM order');
            surf(h,X,Y,sdjln,'FaceColor','interp','EdgeColor','interp');
            zlabel('log norm');
            
            h = pm.nextPlot('jac_deim_diff',...
                'Log norm differences full<->deim','trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,abs(jln-djln));
            zlabel('log norm diff');
            h = pm.nextPlot('jac_st_deim_diff',...
                'Log norm differences full<->simtrans-deim','trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,abs(jln-sdjln));
            zlabel('log norm diff');
            h = pm.nextPlot('rel_diff_full_stdeim',...
                'Relative log norm differences full<->simtrans-deim','trajectory point','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,abs((jln-sdjln)./jln));
            zlabel('rel log norm diff');
            if nargout < 1
                pm.done;
            end
        end
        
        function [ln, aln, s, times] = compareSimTransJac_FullJac(m)
            md = m.Data;
            jtd = md.JacobianTrainingData;
            jstd = md.JacSimTransData;
            
            [~,s] = svd(jstd.VFull,'econ');
            ln = jstd.LogNorms;
            times = jstd.CompTimes;

            aln = ln;
            times_st = times;
            pi = tools.ProcessIndicator('Computing approximated log norms',n*nt);
            for nr = 1:nt
                [x, mu] = md.getTrajectoryNr(nr);    
                for i=1:n
                    J = m.System.f.getStateJacobian(x(:,i),0,mu);
                    t = tic;
                    aln((nr-1)*n+i) = general.Utils.logNorm(v'*J*v);
                    times_st((nr-1)*n+i) = toc(t);
                    pi.step;
                end
            end
            pi.stop;
            times = [times; times_st];
            tmp = mean(times,2);
            fprintf('Mean computation times over %d samples: Full jac log norm %gs, SimTrans log norm %gs\n',...
                n*nt,tmp);
            if nargout == 0
                save JacSimTrans ln aln s times v;
                testing.DEIM.compareSimTransJac_FullJac_plots(ln, aln, s, v);
            end
        end
        
        function pm = compareSimTransJac_FullJac_plots(ln, aln, s, v, pm)
            if nargin < 5
                pm = tools.PlotManager(false,2,2);
            end
            pm.nextPlot('correct ln',sprintf('Full logarithmic norms\nJacobian full dimension %dx%d',...
                size(v,1),size(v,1)));
            plot(ln);
            pm.nextPlot('approx ln',sprintf('Approx logarithmic norms\nSimilarity transformed matrix size %dx%d',...
                size(v,2),size(v,2)));
            plot(aln);
            pm.nextPlot('-','Relative error');
            plot(abs(ln-aln)./ln);
            pm.nextPlot('sing_vals','Singular values of SVD of eigenvector matrix');
            semilogy(diag(s));
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
