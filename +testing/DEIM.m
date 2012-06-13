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
                oldst = deim.ST;
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
            for m = 1:no
                [~, yr, ct] = r.simulate(mu, inidx);
                errs(m,:) = Norm.L2(y-yr);
                relerrs(m,:) = errs(m,:)./Norm.L2(y);
                times(m) = ct;
            end
            d.Order = oldo;
            r.ErrorEstimator.Enabled = olde;
        end
        
        
        %% %%%%%%%%%%%%%%%%%%%% MATRIX DEIM STUFF %%%%%%%%%%%%%%%%%%%%
        function [no,no1,nom,jln,djln,sdjln] = matrix_deim(m, nr)
            deim = m.ErrorEstimator.JacMDEIM;
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
            STFull = m.ErrorEstimator.STFull;
            f = m.System.f;
            
            % Input checks
            if nargin < 2
                st_sizes = 1:size(STFull,2);
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
                    ST = STFull(:,1:st_sizes(snr));
                    t = tic;
                    aln(snr,nr) = general.Utils.logNorm(ST'*J*ST);
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
                m.System.f.XDim,m.System.f.XDim));
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
            
            sv = m.ErrorEstimator.STSingVals;
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
                    jd.setSimilarityTransform(e.STFull(:,1:st_sizes(snr)));
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
