classdef VKOGA
    % VKOGA: Tests for the VKOGA algorithm
    %
    % @author Daniel Wirtz @date 2012-04-04
    %
    % @new{0,6,dw,2012-04-04} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Constant)
        LineWidth = 2;
    end
    
    methods(Static)
        function WH12c(seed)
            %             ap = KerMor.App;
            %             oldv = ap.Verbose;
            %             ap.Verbose = 1;
            if nargin < 1
                seed = 1;
            end
            %dir = '/usr/local/datastore/kermor/VKOGA';
            dir = 'C:/Users/CreaByte/Documents/Uni/DWCAA12/img';
            types = {'jpg'};
            pm = tools.PlotManager;
            pm.LeaveOpen = false;
            pm.UseFileTypeFolders = true;
            pm.NoTitlesOnSave = true;
            
            %             pm.SingleSize = [1224 768];
            %             pm.FilePrefix = 'selection_illus';
            %             %testing.VKOGA.selectCritGraphic(pm);
            %             testing.VKOGA.selectCritGraphic(pm,4554);
            %             pm.done;
            %             pm.savePlots(dir, types, [], true);
            
%             [res, kexp, kexp_oga, a, ao, f, ~] = testing.VKOGA.test_VKOGA_Versions_5dim(1, 1);
%             save tmp res kexp kexp_oga a ao f;
%             return;
%             load tmp;
%             pm.FilePrefix = '5d';
%             pm.SingleSize = [700 800];
%             testing.VKOGA.plotL2Errors(a,ao,pm);
%             pm.SingleSize = [1224 768];
%             testing.VKOGA.plotHerrDecayAndBounds(a,ao,f,pm);
%             pm.done;
%             pm.savePlots(dir, types, [], true);

            runs = 50;
            [res, kexp, kexp_oga, a, ao, f, ~] = testing.VKOGA.test_VKOGA_Versions_5dim(runs, seed);
            save(sprintf('5dim_%druns_%dseed',runs,seed), 'res', 'kexp', 'kexp_oga', 'a', 'ao', 'f');
%             return;
%             load 5dim_runs;
            pm.FilePrefix = '5d_runs';
            pm.SingleSize = [1224 768];
            %pm.SingleSize = [700 800];
            testing.VKOGA.plotVKOGARes(res,pm);
            %testing.VKOGA.plotL2Errors(a,ao,pm);
            %testing.VKOGA.plotHerrDecayAndBounds(a,ao,f,pm);
            pm.done;
            pm.savePlots(dir, types, [], true);
            
%             pm.FilePrefix = '1d';
%             [kexp, kexp_oga, a, ao, f, atd] = testing.VKOGA.test_VKOGA_1D(false, 25);
%             save tmp2 kexp atd a f kexp_oga ao;
%             pm.SingleSize = [1224 768];
%             testing.VKOGA.test_VKOGA_1D_plots(kexp, kexp_oga, a, ao, f, atd, pm);
%             pm.done;
%             pm.savePlots(dir, types, [], true);

            % Dat hier läuft EXTREM schlecht für VKOGA..
%             name = 'Franke3D';%'F7';
%             pm.FilePrefix = sprintf('testfun_%s',name);
%             [kexp, kexpo, a, ao, atd] = testing.VKOGA.test_VKOGA_TestFuns(name, 500, seed);
%             %load(sprintf('vkoga_testfun_%s.mat',name));
%             testing.VKOGA.test_VKOGA_TestFuns_plots(a,ao,pm);
%             pm.LeaveOpen = true;
            
            %             pm = tools.PlotManager(true);
            %             pm.FilePrefix = '5druns';
            %             %load 5dim_VKOGA;
            %             [res, kexp, kexp_OGA, a, ao, f, atd] = testing.VKOGA.test_VKOGA_Versions_5dim(50, 1);
            %             save 5dim_VKOGA res kexp kexp_OGA a ao f atd;
            %             testing.VKOGA.plotVKOGARes(res, pm);
            %             pm.savePlots(dir, 'fig');
            %             pm.savePlots(dir, 'eps', true);
            %

%             % Different gamma values
%             res = testing.VKOGA.test_VKOGA_Versions_diffgamma(30, 1);
%             save diffgamma res;
%             pm.FilePrefix = 'diffgamma';
%             testing.VKOGA.plotVKOGARes(res, pm);
%             pm.done;
%             pm.savePlots(dir, types, [], true);
            
            %             ap.Verbose = oldv;
        end
        
        function pm = selectCritGraphic(pm, seed)
            if nargin < 2
                seed = 6209; %4554, 502197612
                if nargin < 1
                    pm = tools.PlotManager;
                    pm.LeaveOpen = true;
                end
            end
            ms = 16; % marker size
            fine = .01;
            x = -10:fine:10;
            dia = x(end)-x(1);
            r = RandStream('mt19937ar','Seed',seed);
            
            nc = 6;
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            % Use same dist aka space
            f.Kernel.setGammaForDistance(dia/4,1e-5);
            f.Centers.xi = r.rand(1,nc)*dia+x(1);
            f.Ma = r.rand(1,nc)*5-2.5;
            
            fx = f.evaluate(x);
            h = pm.nextPlot(sprintf('seed_%d',seed),'Illustration of VKOGA/WSOGA selection criteria','x');
            plot(h,x,fx,'r-','LineWidth',testing.VKOGA.LineWidth);
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = f.Kernel;
            pos = [-7.1   -5.6   -2.1    1.9    5.9    8.4];
            %pos = [-3 -2];
            c = round((pos-x(1))/fine);
            kexp.Centers.xi = x(c);
            kexp.Ma = (kexp.getKernelMatrix\fx(c)')';
            
            afx = kexp.evaluate(x);
            hold(h,'on');
            plot(h,x,afx,'b-','LineWidth',testing.VKOGA.LineWidth);
            
            % Extension errors
            free = true(size(x));
            free(c) = false;
            xf = x(free);
            Kbig = kexp.getKernelVector(xf)';
            
            % Compute kernel fcn projection to H(m-1)
            A = kexp.getKernelMatrix \ Kbig;
            F = fx(c);
            
            oga_err = abs(fx(free) - F*A);
            %oga_err = abs(fx-afx);
            phinormsq = sqrt(1 - sum(A.*Kbig,1));
            vkoga_err = oga_err ./ phinormsq;
            
            g = [0 .5 0];
            plot(h,xf,oga_err,'Color',g,'LineWidth',testing.VKOGA.LineWidth);
            plot(h,xf,vkoga_err,'--','Color',g,'LineWidth',testing.VKOGA.LineWidth);
            
            off = -3;
            plot(h,xf,phinormsq+off,'m-.','LineWidth',testing.VKOGA.LineWidth);
            
            plot(h,x(c),fx(c),'b.','MarkerSize',ms+8,'LineWidth',testing.VKOGA.LineWidth);
            [om, oidx] = max(oga_err);
            [vm, vidx] = max(vkoga_err);
            plot(h,xf(oidx),om,'rx',xf(vidx),vm,'rx','MarkerSize',ms,'LineWidth',testing.VKOGA.LineWidth);
            
            % Lines at maxima
            plot(h,[xf(oidx) xf(oidx)],[om phinormsq(oidx)+off],'k--');
            plot(h,[xf(vidx) xf(vidx)],[vm phinormsq(vidx)+off],'k--');
            % Zero line
            plot(h,x,0,'k');
            
            lh = legend('f(x)','approx','\langle f_j-f_j^{m-1}, \phi(x,\cdot)\rangle_H',...
                '\langle f_j, \phi^{m-1}_x\rangle_H','||\phi^{m-1}_x||-3');
            set(lh,'Location','Best');
        end
        
        function pm = selectCritGraphic_Diffgamma(pm, seed)
            if nargin < 2
                seed = 6209; %4554, 502197612
                if nargin < 1
                    pm = tools.PlotManager;
                    pm.LeaveOpen = true;
                end
            end
            ms = 16; % marker size
            fine = .01;
            x = -10:fine:10;
            dia = x(end)-x(1);
            r = RandStream('mt19937ar','Seed',seed);
            numgammas = 10;
            gammas = logspace(log10(.01*dia),log10(.1*dia),numgammas);
            truegidx = floor(numgammas/2);
            
            nc = 6;
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            % Use same dist aka space
            f.Kernel.setGammaForDistance(gammas(truegidx),1e-5);
            f.Centers.xi = r.rand(1,nc)*dia+x(1);
            f.Ma = r.rand(1,nc)*5-2.5;
            
            fx = f.evaluate(x);
            h = pm.nextPlot(sprintf('seed_%d',seed),'Illustration of VKOGA/WSOGA selection criteria','x');
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = f.Kernel;
            pos = [-7.1   -5.6   -2.1    1.9    5.9    8.4];
            %pos = [-3 -2];
            c = round((pos-x(1))/fine);
            kexp.Centers.xi = x(c);
            
            % Preps
            free = true(size(x));
            free(c) = false;
            xf = x(free);
            
            %% Compute stuff for "true" approx and choices
            kexp.Kernel.Gamma = gammas(truegidx);
            kexp.Ma = (kexp.getKernelMatrix\fx(c)')';
            tafx = kexp.evaluate(x);
            Kbig = kexp.getKernelVector(xf)';
            A = kexp.getKernelMatrix \ Kbig;
            F = fx(c);
            toga_err = abs(fx(free) - F*A);
            tphinormsq = sqrt(1 - sum(A.*Kbig,1));
            tvkoga_err = toga_err ./ tphinormsq;
            [tom, toidx] = max(toga_err);
            [tvm, tvidx] = max(tvkoga_err);
            
            for k=1:numgammas
                hold(h,'off');
                
                % Full function
                plot(h,x,fx,'r-','LineWidth',testing.VKOGA.LineWidth);
                
                %% Current computation
                kexp.Kernel.Gamma = gammas(k);
                title(h,sprintf('Current Gamma: %g',gammas(k)));
                kexp.Ma = (kexp.getKernelMatrix\fx(c)')';
                afx = kexp.evaluate(x);
                hold(h,'on');
                plot(h,x,afx,'b-','LineWidth',testing.VKOGA.LineWidth);

                % Extension errors
                free = true(size(x));
                free(c) = false;
                xf = x(free);
                Kbig = kexp.getKernelVector(xf)';

                % Compute kernel fcn projection to H(m-1)
                A = kexp.getKernelMatrix \ Kbig;
                F = fx(c);

                oga_err = abs(fx(free) - F*A);
                %oga_err = abs(fx-afx);
                phinormsq = sqrt(1 - sum(A.*Kbig,1));
                vkoga_err = oga_err ./ phinormsq;

                g = [0 .5 0];
                plot(h,xf,oga_err,'Color',g,'LineWidth',testing.VKOGA.LineWidth);
                plot(h,xf,vkoga_err,'--','Color',g,'LineWidth',testing.VKOGA.LineWidth);

                off = -3;
                plot(h,xf,phinormsq+off,'m-.','LineWidth',testing.VKOGA.LineWidth);

                plot(h,x(c),fx(c),'b.','MarkerSize',ms+8,'LineWidth',testing.VKOGA.LineWidth);
                [om, oidx] = max(oga_err);
                [vm, vidx] = max(vkoga_err);
                plot(h,xf(oidx),om,'rx',xf(vidx),vm,'rx','MarkerSize',ms,'LineWidth',testing.VKOGA.LineWidth);

                % Lines at maxima
                plot(h,[xf(oidx) xf(oidx)],[om phinormsq(oidx)+off],'k--');
                plot(h,[xf(vidx) xf(vidx)],[vm phinormsq(vidx)+off],'k--');

                % Zero line
                plot(h,x,0,'k');
                lh = legend('f(x)','approx','\langle f_j-f_j^{m-1}, \phi(x,\cdot)\rangle_H',...
                    '\langle f_j, \phi^{m-1}_x\rangle_H','||\phi^{m-1}_x||-3');
                
                %% "True" space plots
                hold(h,'on');
                plot(h,x,tafx,'-','Color',[.7 .7 1],...
                    'LineWidth',testing.VKOGA.LineWidth);
                tg = [.6 1 .6];
                plot(h,xf,toga_err,'Color',tg,'LineWidth',testing.VKOGA.LineWidth);
                plot(h,xf,tvkoga_err,'--','Color',tg,'LineWidth',testing.VKOGA.LineWidth);
                off = -3;
                plot(h,xf,tphinormsq+off,'-.','Color',[1 .7 1],...
                    'LineWidth',testing.VKOGA.LineWidth);
                plot(h,xf(toidx),tom,'x',xf(tvidx),tvm,'x',...
                    'Color',[1 .7 .7],...
                    'MarkerSize',ms,'LineWidth',testing.VKOGA.LineWidth);

                % Lines at maxima
                %plot(h,[xf(toidx) xf(toidx)],[tom tphinormsq(toidx)+off],'k--');
                %plot(h,[xf(tvidx) xf(tvidx)],[tvm tphinormsq(tvidx)+off],'k--');
                
                set(lh,'Location','Best');
                pause;
            end
        end
        
        function [kexp, kexp_oga, a, ao, f, atd] = test_VKOGA_1D(oga, seed)
            if nargin < 2
                seed = 1;
                if nargin < 1
                    oga = false;
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            dim = 1;
            centers = 12;
            
            x = repmat(-10:.1:10,dim,1);
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.NumGammas = 1;
            a.Dists = sqrt(dim)*20/4;
            a.gameps = .05;
            a.MaxExpansionSize = 100;
            a.UsefScaling = false;
            a.UseOGA = oga;
            a.MaxRelErr = 1e-2;
            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            % Use same dist aka space
            f.Kernel.setGammaForDistance(a.Dists,a.gameps)
            f.Centers.xi = r.rand(dim,centers)*20-10;
            f.Ma = (r.rand(dim,centers)-.5)*150;
            
            fx = f.evaluate(x);
            atd = data.ApproxTrainData(x,[],[]);
            atd.fxi = fx;
            
            vxsize = 1000;
            vx = unique(r.rand(vxsize,dim)*20-10,'rows')';
            a.vxtmuargs{1} = vx;
            a.vfx = f.evaluate(vx);
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            a.computeApproximation(kexp, atd);
            
            ao = a.clone;
            ao.UseOGA = true;
            kexp_oga = kernels.KernelExpansion;
            kexp_oga.Kernel = kernels.GaussKernel;
            ao.computeApproximation(kexp_oga, atd);
        end
        
        function test_VKOGA_1D_plots(kexp, kexp_oga, a, ~, f, atd, pm)
            if nargin < 3
                pm = tools.PlotManager(false, 1, 2);
                pm.LeaveOpen = true;    
            end
            
%             x = repmat(-15:.01:15,1,1);
%             fx = f.evaluate(x);
            x = atd.xi.toMemoryMatrix;
            fx = atd.fxi.toMemoryMatrix;
            %afx = kexp.evaluate(x);
            %err2 = Norm.L2((fx-afx)');
            h = pm.nextPlot('vkoga_approx',sprintf('Max L^2-error on training data: %g. Centers: VKOGA %d, WSOGA2 %d',...
                a.MaxRelErr,size(kexp.Centers.xi,2),size(kexp_oga.Centers.xi,2)));
            plot(h,x',fx','k--',...
                f.Centers.xi, f.evaluate(f.Centers.xi), 'kx',...
                kexp.Centers.xi,f.evaluate(kexp.Centers.xi),'ro',...
                kexp_oga.Centers.xi,f.evaluate(kexp_oga.Centers.xi),'bd',...
                'MarkerSize',12,'LineWidth',testing.VKOGA.LineWidth);
            legend(h,'f(x)','f(x) Centers','VKOGA Centers','WSOGA2 Centers');
            
            h = pm.nextPlot('vkoga_bounds');
            n = a.expsizes(1);
            M = f.MBnd;
            bound = M^2 ./ (1+(1:n));
            semilogy(h,1:n,f.NativeNorm^2 - a.HerrDecay(1:n),'r',1:n, bound,'b',...
                'LineWidth',testing.VKOGA.LineWidth);%,1:n, bound2,'b--');
            legend(h,'H-approximation error at step m','Theoretical upper bound');
        end
        
        function [res, kexp, kexp_OGA, a, ao, f, atd] = test_VKOGA_Versions_diffgamma(runs, seed)
            % res structure: 8 x runs double with
            % 1-2: absolute errors on training set for komp,oga
            % 3-4: absolute errors on validation set for komp,oga
            % 5-6: expansion sizes for komp,oga when achieving MaxRelErr
            % 7-8: relative errors on validation set for komp,oga
            if nargin < 2
                seed = 1;
                if nargin < 1
                    runs = 100;
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            rf = RandStream('mt19937ar','Seed',seed^2 + 1);
            centers = 20;
            xrange = 10; xoff = -5;
            dim = 5;
            atdsize = max(dim*500,1000);
            vxsize = max(dim*200,500);
            
            x = r.rand(dim,atdsize)*xrange+xoff;
            if dim == 1
                x = sort(x);
            end
            xu = unique(x','rows')';
            if size(xu,2) ~= size(x)
                x = xu;
                warning('crap:id','Created nonunique centers!');
            end
            dia = sqrt(dim)*xrange;
            %             % Get realistic training data
            %             m = models.pcd.PCDModel;
            %             mu = 0.0081;%m.System.getRandomParam;
            %             m.T = 20;
            %             [~,x] = m.simulate(mu);
            %             x = x(1:2,:);
            %             dim = size(x,1);
            %             [m, M] = general.Utils.getBoundingBox(x);
            %             xrange = max(M-m); xoff = min(m);
            %             dia = Norm.L2(M-m);
            
            dists = linspace(dia/15,dia*2,runs);
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,200);
            a.UsefScaling = false;
            a.MaxRelErr = 1e-3;
            a.PhiNormMin = sqrt(eps);%1e-4;%
            
            vx = unique(r.rand(vxsize,dim)*xrange+xoff,'rows')';
            a.vxtmuargs{1} = vx;
            
            ao = a.clone;
            ao.UseOGA = true;
            
            % Setup test functions
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            res = struct;
            res.a = a;
            res.dists = dists;
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                % Choose same but run-dependent gamma values
                a.Dists = dists(run);
                ao.Dists = dists(run);
                f.Kernel.setGammaForDistance(dists(run),a.gameps);
                res.gammas(run) = f.Kernel.Gamma;
                
                f.Centers.xi = rf.rand(dim,centers)*xrange + xoff;
                f.Ma = rf.rand(dim,centers)*15;%-2.5;
                
                atd.fxi = f.evaluate(x);
                
                % Set validation f-evaluation
                a.f = f;
                ao.f = f;
                a.vfx = f.evaluate(vx);
                ao.vfx = a.vfx;
                
                % Original greedy version
                ao.computeApproximation(kexp, atd);
                kexp_OGA = kexp.clone;
                % Normalized version
                a.computeApproximation(kexp, atd);
                
                % Store error information at points of same expansion
                % i.e. subspace size!
                s = min(size(kexp.Ma,2),size(kexp_OGA.Ma,2));
                res.terro(:,run) = ao.err(:,s);
                res.verro(:,run) = ao.verr(:,s);
                res.Herro(:,run) = f.NativeNorm^2 - ao.HerrDecay(:,s);
                res.vrelerro(:,run) = ao.vrelerr(:,s);
                res.cnumo(run) = size(kexp_OGA.Ma,2);
                res.ogabnd(run) = f.MBnd^2 * dim / (1 + s/dim);
                
                res.terr(:,run) = a.err(:,s);
                res.verr(:,run) = a.verr(:,s);
                res.Herr(:,run) = f.NativeNorm^2 - a.HerrDecay(:,s);
                res.vrelerr(:,run) = a.vrelerr(:,s);
                res.cnum(run) = size(kexp.Ma,2);
                res.vkogabnd(run) = f.MBnd^2 * a.VKOGABound(s);
                
                
                pi.step(run);
                
                if runs == 1 || true
                    testing.VKOGA.plotStatistics(a,ao,f);
                    if runs > 1
                        %pause;
                    end
                end
            end
            pi.stop;
            
            if nargout == 0
                testing.VKOGA.plotKompRes(res);
            end
        end
        
        function [res, kexp, kexp_OGA, a, ao, f, atd] = test_VKOGA_Versions_tube(dim, runs, seed)
            if nargin < 3
                seed = 1;
                if nargin < 2
                    runs = 100;
                    if nargin < 1
                        dim = 5;
                    end
                end
            end
            rf = RandStream('mt19937ar','Seed',2*seed);
            centers = 20;
            xrange = 30;
            atdsize = 2500;
            vxsize = 200;
            
            spread = .15/3;
            x = general.Utils.getTube(dim, atdsize+vxsize+centers*runs, xrange, spread, 3*seed);
            % Take away the centers
            xi = x(:,1:centers*runs);
            x(:,1:centers*runs) = [];
            dia = spread*xrange/5;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            a.Dists = dia;
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,400);
            a.UsefScaling = false;
            a.MaxRelErr = 1e-3;
            if dim == 1
                a.PhiNormMin = .0001;
            else
                a.PhiNormMin = sqrt(eps);
            end
            
            %vx = unique(r.rand(vxsize,dim)*xrange+xoff,'rows')';
            vx = x(:,1:vxsize);
            a.vxtmuargs{1} = vx;
            x(:,1:vxsize) = [];
            
            ao = a.clone;
            ao.UseOGA = true;
            
            % Setup test functions
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(dia,a.gameps)
            
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            res = struct;
            res.a = a;
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                f.Centers.xi = xi(:,1+(run-1)*centers:run*centers);
                f.Ma = rf.rand(dim,centers)*15;%-2.5;
                atd.fxi = f.evaluate(x);
                
                % Set validation f-evaluation
                a.f = f;
                ao.f = f;
                a.vfx = f.evaluate(vx);
                ao.vfx = a.vfx;
                
                % Original greedy version
                ao.computeApproximation(kexp, atd);
                kexp_OGA = kexp.clone;
                % Normalized version
                a.computeApproximation(kexp, atd);
                
                % Store error information at points of same expansion
                % i.e. subspace size!
                s = min(size(kexp.Ma,2),size(kexp_OGA.Ma,2));
                res.terro(:,run) = ao.err(:,s);
                res.verro(:,run) = ao.verr(:,s);
                res.Herro(:,run) = f.NativeNorm^2 - ao.HerrDecay(:,s);
                res.vrelerro(:,run) = ao.vrelerr(:,s);
                res.cnumo(run) = size(kexp_OGA.Ma,2);
                res.ogabnd(run) = f.MBnd^2 * dim / (1 + s/dim);
                
                res.terr(:,run) = a.err(:,s);
                res.verr(:,run) = a.verr(:,s);
                res.Herr(:,run) = f.NativeNorm^2 - a.HerrDecay(:,s);
                res.vrelerr(:,run) = a.vrelerr(:,s);
                res.cnum(run) = size(kexp.Ma,2);
                res.vkogabnd(run) = f.MBnd^2 * a.VKOGABound(s);
                
                pi.step;
                
                if runs == 1
                    testing.VKOGA.plotStatistics(a,ao,f);
                    if runs > 1
                        %pause;
                    end
                end
            end
            pi.stop;
            
            if nargout == 0 && runs > 1
                testing.VKOGA.plotVKOGARes(res);
            end
        end
        
        function [res, kexp, kexp_OGA, a, ao, f, atd] = test_tube_fcomp_diffcenter(dim, runs, seed)
            if nargin < 3
                seed = 1;
                if nargin < 2
                    runs = 100;
                    if nargin < 1
                        dim = 5;
                    end
                end
            end
            rf = RandStream('mt19937ar','Seed',2*seed);
            compcenters = 10;
            centers = dim*compcenters;
            xrange = 10;
            atdsize = 2500;
            vxsize = 2000;
            
            spread = .15;
            x = general.Utils.getTube(dim, atdsize+vxsize+centers*runs, xrange, spread, 3*seed);
            % Take away the centers
            xi = x(:,1:centers*runs);
            x(:,1:centers*runs) = [];
            dia = spread*xrange/2;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            a.Dists = dia;
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,800);
            a.UsefScaling = false;
            a.MaxRelErr = 1e-2;
            a.PhiNormMin = sqrt(eps);
            
            %vx = unique(r.rand(vxsize,dim)*xrange+xoff,'rows')';
            vx = x(:,1:vxsize);
            a.vxtmuargs{1} = vx;
            x(:,1:vxsize) = [];
            
            ao = a.clone;
            ao.UseOGA = true;
            
            % Setup test functions
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(dia,a.gameps)
            
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            res = struct;
            res.a = a;
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                f.Centers.xi = xi(:,1+(run-1)*centers:run*centers);
                f.Ma = zeros(dim,centers);
                % set 20 effective centers for each component
                for di = 1:dim
                    f.Ma(di,(1+(di-1)*compcenters):di*compcenters) = rf.rand(1,compcenters)*15;
                end
                atd.fxi = f.evaluate(x);
                
                % Set validation f-evaluation
                a.f = f;
                ao.f = f;
                a.vfx = f.evaluate(vx);
                ao.vfx = a.vfx;
                
                % Original greedy version
                ao.computeApproximation(kexp, atd);
                kexp_OGA = kexp.clone;
                % Normalized version
                a.computeApproximation(kexp, atd);
                
                % Store error information at points of same expansion
                % i.e. subspace size!
                s = min(size(kexp.Ma,2),size(kexp_OGA.Ma,2));
                res.terro(:,run) = ao.err(:,s);
                res.verro(:,run) = ao.verr(:,s);
                res.Herro(:,run) = f.NativeNorm^2 - ao.HerrDecay(:,s);
                res.vrelerro(:,run) = ao.vrelerr(:,s);
                res.cnumo(run) = size(kexp_OGA.Ma,2);
                res.ogabnd(run) = f.MBnd^2 * dim / (1 + s/dim);
                
                res.terr(:,run) = a.err(:,s);
                res.verr(:,run) = a.verr(:,s);
                res.Herr(:,run) = f.NativeNorm^2 - a.HerrDecay(:,s);
                res.vrelerr(:,run) = a.vrelerr(:,s);
                res.cnum(run) = size(kexp.Ma,2);
                res.vkogabnd(run) = f.MBnd^2 * a.VKOGABound(s);
                
                pi.step;
                
                if runs == 1 || true
                    testing.VKOGA.plotStatistics(a,ao,f);
                    if runs > 1
                        %pause;
                    end
                end
            end
            pi.stop;
            
            if nargout == 0
                testing.VKOGA.plotKompRes(res);
            end
        end
        
        function [res, kexp, kexp_OGA, a, ao, f, atd] = test_VKOGA_Versions_5dim(runs, seed)
            % res structure: 8 x runs double with
            % 1-2: absolute errors on training set for komp,oga
            % 3-4: absolute errors on validation set for komp,oga
            % 5-6: expansion sizes for komp,oga when achieving MaxRelErr
            % 7-8: relative errors on validation set for komp,oga
            if nargin < 2
                seed = 1;
                if nargin < 1
                    runs = 100;
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            rf = RandStream('mt19937ar','Seed',seed^2 + 1);
            centers = 20;
            xrange = 10; xoff = -5;
            dim = 5;
            atdsize = max(dim*500,1000);
            vxsize = max(dim*200,500);
            
            x = r.rand(dim,atdsize)*xrange+xoff;
            if dim == 1
                x = sort(x);
            end
            xu = unique(x','rows')';
            if size(xu,2) ~= size(x)
                x = xu;
                warning('crap:id','Created nonunique centers!');
            end
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            if runs == 1
                dia = sqrt(dim)*xrange / 2.5;
                a.gameps = .4;
                a.MaxRelErr = 1e-3;
            else
                % "Old" settings for multiple runs
                a.gameps = .6;
                a.MaxRelErr = 1e-4;
                dia = sqrt(dim*xrange);
            end
            a.Dists = dia;
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,200);
            a.UsefScaling = false;
            a.PhiNormMin = sqrt(eps);
            
            vx = unique(r.rand(vxsize,dim)*xrange+xoff,'rows')';
            a.vxtmuargs{1} = vx;
            
            ao = a.clone;
            ao.UseOGA = true;
            
            % Setup test functions
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(dia,a.gameps)
            
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            res = struct;
            res.a = a;
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                f.Centers.xi = rf.rand(dim,centers)*xrange + xoff;
                f.Ma = rf.rand(dim,centers)*15;%-2.5;
                atd.fxi = f.evaluate(x);
                
                % Set validation f-evaluation
                a.f = f;
                ao.f = f;
                a.vfx = f.evaluate(vx);
                ao.vfx = a.vfx;
                
                % Original greedy version
                ao.computeApproximation(kexp, atd);
                kexp_OGA = kexp.clone;
                % Normalized version
                a.computeApproximation(kexp, atd);
                
                % Store error information at points of same expansion
                % i.e. subspace size!
                s = min(size(kexp.Ma,2),size(kexp_OGA.Ma,2));
                res.terro(:,run) = ao.err(:,s);
                res.verro(:,run) = ao.verr(:,s);
                res.Herro(:,run) = f.NativeNorm^2 - ao.HerrDecay(:,s);
                res.vrelerro(:,run) = ao.vrelerr(:,s);
                res.cnumo(run) = size(kexp_OGA.Ma,2);
                res.ogabnd(run) = f.MBnd^2 * dim / (1 + s/dim);
                
                res.terr(:,run) = a.err(:,s);
                res.verr(:,run) = a.verr(:,s);
                res.Herr(:,run) = f.NativeNorm^2 - a.HerrDecay(:,s);
                res.vrelerr(:,run) = a.vrelerr(:,s);
                res.cnum(run) = size(kexp.Ma,2);
                res.vkogabnd(run) = f.MBnd^2 * a.VKOGABound(s);
                
                pi.step;
            end
            pi.stop;
        end
        
        function [res, kexp, kexp1, a, a1, f, atd] = test_VKOGA_Versions_5dim_2step(runs, seed)
            % res structure: 8 x runs double with
            % 1-2: absolute errors on training set for komp,oga
            % 3-4: absolute errors on validation set for komp,oga
            % 5-6: expansion sizes for komp,oga when achieving MaxRelErr
            % 7-8: relative errors on validation set for komp,oga
            if nargin < 2
                seed = 1;
                if nargin < 1
                    runs = 100;
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            rf = RandStream('mt19937ar','Seed',seed^2 + 1);
            centers = 10;
            xrange = 10; xoff = -5;
            dim = 5;
            atdsize = 300;
            vxsize = 500;
            
            x = r.rand(dim,atdsize)*xrange+xoff;
            if dim == 1
                x = sort(x);
            end
            xu = unique(x','rows')';
            if size(xu,2) ~= size(x)
                x = xu;
                warning('crap:id','Created nonunique centers!');
            end
            %             x = repmat(xoff+(0:xrange/10:xrange),dim,1);
            %[m,M] = general.Utils.getBoundingBox(x);
            %             dia = Norm.L2(M-m);
            %             dia = sqrt(dim)*xrange;
            dist = sqrt(dim*xrange);
            
            a = approx.algorithms.VectorialKernelOMP2Step;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            a.Dists = dist;
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,200);
            a.UsefScaling = false;
            a.MaxRelErr = 1e-3;
            a.PhiNormMin = sqrt(eps);
            
            a1 = approx.algorithms.VectorialKernelOMP;
            a1.CoeffComp = general.interpolation.KernelInterpol;
            a1.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            a1.Dists = a.Dists;
            a1.NumGammas = a.NumGammas;
            a1.MaxExpansionSize = a.MaxExpansionSize;
            a1.UsefScaling = a.UsefScaling;
            a1.MaxRelErr = a.MaxRelErr;
            a1.PhiNormMin = a.PhiNormMin;
            
            vx = unique(r.rand(vxsize,dim)*xrange+xoff,'rows')';
            a.vxtmuargs{1} = vx;
            a1.vxtmuargs{1} = vx;
            
            % Setup test functions
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(dist,a.gameps)
            
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            res = struct;
            res.a = a;
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                %f.Centers.xi = x(:,1:centers);
                f.Centers.xi = rf.rand(dim,centers)*xrange + xoff;
                f.Ma = rf.rand(dim,centers)*15;%-2.5;
                atd.fxi = f.evaluate(x);
                
                % Set validation f-evaluation
                a.f = f;
                a1.f = f;
                a.vfx = f.evaluate(vx);
                a1.vfx = a.vfx;
                
                % Original VKOGA
                a.computeApproximation(kexp, atd);
                kexp1 = kexp.clone;
                % 2-Step
                a1.computeApproximation(kexp1, atd);
                
                % Store error information at points of same expansion
                % i.e. subspace size!
                s = min(size(kexp.Ma,2),size(kexp1.Ma,2));
                %s = size(kexp.Ma,2);
                res.terro(:,run) = a1.err(:,s);
                res.verro(:,run) = a1.verr(:,s);
                res.Herro(:,run) = f.NativeNorm^2 - a1.HerrDecay(:,s);
                res.vrelerro(:,run) = a1.vrelerr(:,s);
                res.cnumo(run) = size(kexp1.Ma,2);
                res.ogabnd(run) = f.MBnd^2 * dim / (1 + s/dim);
                
                res.terr(:,run) = a.err(:,s);
                res.verr(:,run) = a.verr(:,s);
                res.Herr(:,run) = f.NativeNorm^2 - a.HerrDecay(:,s);
                res.vrelerr(:,run) = a.vrelerr(:,s);
                res.cnum(run) = size(kexp.Ma,2);
                res.vkogabnd(run) = f.MBnd^2 * a.VKOGABound(s);
                
                pi.step(run);
                
                if runs == 1
                    %testing.VKOGA.plotStatistics(a,a1,f);
                    if runs > 1
                        %pause;
                    end
                end
            end
            pi.stop;
            
            if nargout == 0 && runs > 1
                testing.VKOGA.plotKompRes(res);
            end
        end
        
        function test_Gain
            % Tests if the gain computed by the scalar product is equal to the gain directly
            % computed by comparison of subsequent interpolations after corresponding subspace
            % extension.
            dim = 1;
            centers = 3;
            atdsize = 50;
            r = RandStream('mt19937ar','Seed',1);
            
            x = r.rand(dim,atdsize)*10-5;
            %if dim > 1
            x = unique(x','rows')';
            %end
            %             [m,M] = general.Utils.getBoundingBox(x);
            %             dia = sqrt(sum((M-m).^2));
            %             dis = sqrt(dia);
            dis = 2;
            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(dis,sqrt(eps))
            f.Centers.xi = r.rand(dim,centers)*20-10;
            f.Ma = r.rand(dim,centers)*15;%-2.5;
            nativenorm = f.NativeNorm
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            fxi = f.evaluate(x);
            
            i = general.interpolation.KernelInterpol;
            
            for s = 1:atdsize-1
                used = 1:s;
                free = (s+1):atdsize;
                kexp.Centers.xi = x(:,used);
                K = kexp.getKernelMatrix;
                i.init(data.MemoryKernelMatrix(K));
                
                Ma1 = i.interpolate(fxi(:,used));
                
                kexp.Ma = Ma1';
                fx = kexp.evaluate(x);
                
                base = -2*fxi(:,used)*Ma1;
                base = base + Ma1'*(K*Ma1);
                
                % Compute projections
                Kbig = kexp.getKernelVector(x(:,free))';
                A = i.interpolate(Kbig');
                F = fxi(:,used);
                
                fDotPhiSqAll = (fxi(:,free) - F*A).^2;
                phinormsq = 1 - sum(A.*Kbig,1);
                
                % Erweiterung
                extparts = zeros(size(free));
                cnt=1;
                for ext=free
                    usede = [used ext];
                    kexp.Centers.xi = x(:,usede);
                    K = kexp.getKernelMatrix;
                    i.init(data.MemoryKernelMatrix(K));
                    
                    Ma2 = i.interpolate(fxi(:,usede));
                    
                    hlp = -2*fxi(:,usede)*Ma2;
                    hlp = hlp + Ma2'*(K*Ma2);
                    extparts(cnt) = hlp;
                    cnt=cnt+1;
                end
                
                dif = base - extparts;
                fDotPhiSqAll = fDotPhiSqAll ./ phinormsq;
                figure(1);
                subplot(2,3,1);
                plot(fDotPhiSqAll); title('fDotPhiSqAll'); axis tight;
                subplot(2,3,2);
                plot(dif - fDotPhiSqAll); title('difference');  axis tight;
                subplot(2,3,3);
                semilogy(phinormsq); title('phinormsq');  axis tight;
                subplot(2,3,4);
                plot(nativenorm^2 + extparts); title('nativenorm - projection parts of ||f-I_X f||^2');  axis tight;
                subplot(2,3,5);
                plot(x,f.evaluate(x),'r',x,fx,'b',x(used),fx(used),'b.','MarkerSize',15);
                pause;
            end
        end
        
        function plotL2Errors(a, ao, pm)
            n1 = a.expsizes;
            n2 = ao.expsizes;
            rows = 1;
            if ~isempty(a.vxtmuargs)
                rows = 2;
            end
            if nargin < 3
                pm = tools.PlotManager(false, rows, 3);
            end
            
            h = pm.nextPlot('abs_l2_err_train','Absolute L^2 errors on training set',...
                'expansion size N','absolute error');
            semilogy(h,1:n1,a.err(1:n1),'r',1:n2,ao.err(1:n2),'b','LineWidth',testing.VKOGA.LineWidth);
            legend(h,'VKOGA','WSOGA2');
            
            h = pm.nextPlot('rel_l2_err_train','Relative L^2 errors on training set',...
                'expansion size N','relative error');
            semilogy(h,1:n1,a.relerr(1:n1),'r',1:n2,ao.relerr(1:n2),'b','LineWidth',testing.VKOGA.LineWidth);
            legend(h,'VKOGA','WSOGA2');
            
            if ~isempty(a.vxtmuargs)
                %vsiz = size(a.vfx,2);
                h = pm.nextPlot('abs_l2_err_val','Absolute L^2 errors on validation set',...
                    'expansion size N','absolute errors');
                semilogy(h,1:n1,a.verr(1:n1),'r',1:n2,ao.verr(1:n2),'b','LineWidth',testing.VKOGA.LineWidth);
                legend('VKOGA','WSOGA2');
                
                h = pm.nextPlot('rel_l2_err_val','Relative L^2 errors on validation set',...
                    'expansion size N','relative errors');
                semilogy(h,1:n1,a.vrelerr(1:n1),'r',1:n2,ao.vrelerr(1:n2),'b','LineWidth',testing.VKOGA.LineWidth);
                legend('VKOGA','WSOGA2');
            end
        end
        
        function plotHerrDecayAndBounds(a, ao, f, pm)
            n1 = a.expsizes;
            n2 = ao.expsizes;
            if nargin < 4
                pm = tools.PlotManager;
            end
            
            M = f.MBnd;
            h = pm.nextPlot('H_err_decay_bounds',sprintf('H-norm decay plus upper convergence bound with M=%3.2f', M),...
                'expansion size N');
            fnsq = f.NativeNorm^2;
            semilogy(h,0:n1,fnsq-[0 a.HerrDecay(1:n1)],'r',0:n2,fnsq-[0 ao.HerrDecay(1:n2)],'b','LineWidth',testing.VKOGA.LineWidth);
            hold on;
            dim = size(f.Ma,1);
            bound_sq = (dim*M^2) ./ (1 + (0:(max(n1,n2))) / dim);
            semilogy(h,0:max(n1,n2),bound_sq(1:max(n1,n2)+1),'b--','LineWidth',testing.VKOGA.LineWidth);
            bound_cm = M^2 * [dim a.VKOGABound];
            semilogy(h,0:n1,bound_cm(1:n1+1),'r--','LineWidth',testing.VKOGA.LineWidth);
            legend('VKOGA','WSOGA2','Upper bnd','VKOGA a-post. bnd');
            hold off;
        end
        
        function plotVKOGARes(res, pm)
            if nargin < 2
                pm = tools.PlotManager(false, 2, 3);
                pm.FilePrefix = 'plotVKOGARes';
            end
            lw = testing.VKOGA.LineWidth;
            runs = length(res.terr);
            runx = 1:runs;
            VKOGA_b = res.verr <= res.verro;
            h = pm.nextPlot('errcomp',sprintf('Maximum L^2-errors on training/validation set\nVKOGA <= WSOGA2 on validation set: %2.2f%%',...
                100*sum(VKOGA_b)/runs),'test run','absolute L2 errors');
            plot(h, runx, [res.terr; res.terro; res.verr; res.verro],'LineWidth',lw);
            legend('VKOGA','WSOGA2','VKOGA (val)','WSOGA2 (val)');
            
            VKOGA_b = res.vrelerr <= res.vrelerro;
            pm.nextPlot('relerrcomp',sprintf('Maximum relative L^2-errors on validation set\nVKOGA <= WSOGA2: %2.2f%%',...
                100*sum(VKOGA_b)/runs),'test run','relative L2 errors');
            plot(runx, [res.vrelerr; res.vrelerro],'LineWidth',lw);
            legend('VKOGA (val)','WSOGA2 (val)');
            
            VKOGA_b = res.cnum <= res.cnumo;
            h = pm.nextPlot('expsizes',sprintf('Expansion sizes at MaxRelErr=%e\nVKOGA better than WSOGA2: %2.2f%%',...
                res.a.MaxRelErr,100*sum(VKOGA_b)/runs),'test run','expansion sizes');
            plot(h,runx, [res.cnum; res.cnumo],'LineWidth',lw); 
%             hold on;
%             plot(h,runx(VKOGA_b),res.cnum(VKOGA_b),'g.','MarkerSize',20);
%             plot(h,runx(~VKOGA_b),res.cnumo(~VKOGA_b),'r.','MarkerSize',20);
            legend('VKOGA','WSOGA2');
            
            VKOGA_b = res.Herr <= res.Herro;
            h = pm.nextPlot('Herr_and_bounds',sprintf('H-errors and bounds at MaxRelErr=%e\nVKOGA better than WSOGA2: %2.2f%%',...
                res.a.MaxRelErr,100*sum(VKOGA_b)/runs),'test run','H-errors');
            semilogy(h,runx, res.ogabnd, 'b--', runx, res.vkogabnd, 'r--');
            hold on;
            semilogy(h,runx, res.Herro, 'b', runx, res.Herr, 'r');
            semilogy(h,runx(VKOGA_b),res.Herr(VKOGA_b),'g.','MarkerSize',20);
            semilogy(h,runx(~VKOGA_b),res.Herro(~VKOGA_b),'r.','MarkerSize',20);
            legend('VOGA bound (a-pri)','VKOGA bound (a-post)','H-err WSOGA2','H-err VKOGA');
            
            h = pm.nextPlot('errbnd_ratios','WSOGA2 to VKOGA ratios');
            hlp = res.ogabnd ./ res.vkogabnd;
            plot(h,runx, hlp, 'r', runx, repmat(mean(hlp),1,runs), 'r--');
            hold on;
            hlp = res.Herro ./ res.Herr;
            plot(h,runx, hlp, 'b', runx, repmat(mean(hlp),1,runs), 'b--');
            legend('WSOGA2 to VKOGA bound ratio','average','WSOGA2 to VKOGA Herr ratio','average');
            
            if isfield(res,'gammas')
                pm.nextPlot('gammas');
                plot(1:length(res.gammas),res.gammas);
                legend('gamma values');
                xlabel('run');
                ylabel('gamma value');
            end
        end
        
        function plotErrorComp(vkoga, oga, pm)
            if nargin < 3
                pm = tools.PlotManager(false,1,2);
            end
            h1 = pm.nextPlot('abs_err','Absolute errors','expansion size','error');
            h2 = pm.nextPlot('rel_err','Relative errors','expansion size','error');
            pm.done;
            for k=1:size(vkoga.err,1)
                semilogy(h1,[vkoga.err(k,:); oga.err(k,:)]');
                semilogy(h2,[vkoga.relerr(k,:); oga.relerr(k,:)]');
                title(h1,sprintf('Gammacomp-dist: %g',vkoga.Dists(k)));
                legend(h1,'VKOGA','WSOGA2');
                pause;
            end
        end
        
        function [kexp, kexpo, a, ao, atd] = test_VKOGA_TestFuns(name, atdsize, seed)
            if nargin < 3
                seed = 1;
                if nargin < 2
                    atdsize = 4000;
                    if nargin < 1
                        name = 'F3';
                    end
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            [fun, dom] = eval(sprintf('testing.TestFunctions.%s',name));
            dim = size(dom,1);
            vxsize = atdsize/2;
            
            x = repmat(dom(:,1),1,atdsize+vxsize) + r.rand(dim,atdsize+vxsize).*repmat(dom(:,2)-dom(:,1),1,atdsize+vxsize);
            if dim == 1
                x = sort(x);
            end
            
            a = eval(['testing.VKOGA.getAlg' name]);
            a.f = fun;
            
            %             vx = x(:,1:vxsize);
            %             x(:,1:vxsize) = [];
            %             a.vxtmuargs{1} = vx;
            
            ao = a.clone;
            ao.UseOGA = true;
            ao.f = fun;
            
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            atd.fxi = fun(x);
            % Original greedy version
            ao.computeApproximation(kexp, atd);
            kexpo = kexp.clone;
            % Normalized version
            a.computeApproximation(kexp, atd);
            
            save(sprintf('vkoga_testfun_%s',name),'a','ao','kexp','kexpo','atd');
        end
        
        function test_VKOGA_TestFuns_plots(a, ao, pm)
            pm.nextPlot('VKOGA_abserr','VKOGA absolute error','expansion size','absolute error');
            semilogy(a.err');
            pm.nextPlot('WSOGA2_abserr','WSOGA2 absolute error','expansion size','absolute error');
            semilogy(ao.err');
            pm.nextPlot('VKOGA_relerr','VKOGA relative error','expansion size','relative error');
            semilogy(a.relerr');
            hold on;
            plot(1:a.MaxExpansionSize,a.MaxRelErr,'k-');
            pm.nextPlot('WSOGA2_relerr','WSOGA2 relative error','expansion size','relative error');
            semilogy(ao.relerr');
            hold on;
            plot(1:ao.MaxExpansionSize,ao.MaxRelErr,'k-');
            
            pm.nextPlot('expsizes','Expansion sizes','gamma value','size');
            plot([a.expsizes; ao.expsizes]');
            legend('VKOGA','WSOGA2');
            
            pm.nextPlot('VKOGA_herrdecay','VKOGA H-error decay','iteration nr','decay');
            semilogy(a.HerrDecay');
            pm.nextPlot('WSOGA2_herrdecay','WSOGA2 H-error decay','iteration nr','decay');
            semilogy(ao.HerrDecay');
        end
    end
    
    methods(Static, Access=private)
        function a = getAlgF7
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            %a.Dists = sqrt(dia);
            a.gameps = 1e-3;
            a.NumGammas = 5;
            a.MaxExpansionSize = 200;
            a.UsefScaling = false;
            a.MaxRelErr = 1e-2;
            a.PhiNormMin = sqrt(eps);%1e-3;
            a.MinGFactor = .1;
            a.MaxGFactor = 2;
        end
        
        function a = getAlgFranke3D
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            %a.Dists = sqrt(dia);
            a.gameps = 1e-3;
            a.NumGammas = 10;
            a.MaxExpansionSize = 250;
            a.UsefScaling = false;
            a.MaxRelErr = 1e-2;
            a.PhiNormMin = sqrt(eps);%1e-3;
            a.MinGFactor = .1;
            a.MaxGFactor = 1;
        end
        
        function a = getAlgF3
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            %a.Dists = sqrt(dia);
            a.gameps = 1e-3;
            a.NumGammas = 20;
            a.MaxExpansionSize = 200;
            a.UsefScaling = false;
            a.MaxRelErr = 1e-4;
            a.PhiNormMin = 1e-3;%sqrt(eps);
            a.MinGFactor = .1;
            a.MaxGFactor = 1;
        end
    end
    
end