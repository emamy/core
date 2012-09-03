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
    
    methods(Static)
        function WH12c
%             ap = KerMor.App;
%             oldv = ap.Verbose;
%             ap.Verbose = 1;
            
            dir = '/usr/local/datastore/kermor/VKOGA';
            types = {'jpg'};
            pm = tools.PlotManager;
            pm.SingleSize = [1224 768];
            pm.LeaveOpen = false;
            pm.UseFileTypeFolders = true;
            pm.NoTitlesOnSave = true;
            
            pm.FilePrefix = 'selection_illus';
            testing.VKOGA.selectCritGraphic(pm);
            %testing.VKOGA.selectCritGraphic(pm,4554);
            pm.done;
            pm.savePlots(dir, types, [], true);
            
            %dir = 'C:\Users\CreaByte\Documents\Uni\VKOGA\img\test';
%             [~, ~, ~, a, ao, f, ~] = testing.VKOGA.test_VKOGA_Versions_5dim(1, 1);
%             pm = tools.PlotManager(true);
%             pm.FilePrefix = '5d';
%             testing.VKOGA.plotStatistics(a,ao,f,pm);
%             pm.savePlots(dir, 'fig');
%             pm.savePlots(dir, 'eps', true);
%             
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
%             pm = tools.PlotManager(true);
%             pm.FilePrefix = 'diffgamma';
%             testing.VKOGA.plotVKOGARes(res, pm);
%             pm.savePlots(dir, 'fig');
%             pm.savePlots(dir, 'eps', true);
            
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
            lw = 2;
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
            plot(h,x,fx,'r-','LineWidth',lw);
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = f.Kernel;
            pos = [-7.1   -5.6   -2.1    1.9    5.9    8.4];
            %pos = [-3 -2];
            c = round((pos-x(1))/fine);
            kexp.Centers.xi = x(c);
            kexp.Ma = (kexp.getKernelMatrix\fx(c)')';
            
            afx = kexp.evaluate(x);
            hold(h,'on');
            plot(h,x,afx,'b-','LineWidth',lw);
            
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
            plot(h,xf,oga_err,'Color',g,'LineWidth',lw);
            plot(h,xf,vkoga_err,'--','Color',g,'LineWidth',lw);
            
            off = -3;
            plot(h,xf,phinormsq+off,'m-.','LineWidth',lw);
            
            
            plot(h,x(c),fx(c),'b.','MarkerSize',ms+8,'LineWidth',lw);
            [om, oidx] = max(oga_err);
            [vm, vidx] = max(vkoga_err);
            plot(h,xf(oidx),om,'rx',xf(vidx),vm,'rx','MarkerSize',ms,'LineWidth',lw);
            
            % Lines at maxima
            plot(h,[xf(oidx) xf(oidx)],[om phinormsq(oidx)+off],'k--');
            plot(h,[xf(vidx) xf(vidx)],[vm phinormsq(vidx)+off],'k--');
            % Zero line
            plot(h,x,0,'k');
            
            lh = legend('f(x)','approx','\langle f_j-f_j^{m-1}, \phi(x,\cdot)\rangle_H',...
                '\langle f_j, \phi^{m-1}_x\rangle_H','||\phi^{m-1}_x||-3');
            set(lh,'Location','Best');
        end
        
        function [kexp, atd, a, f] = test_VKOGA(oga, seed, pm)
            if nargin < 3
                pm = tools.PlotManager(false, 1, 2);
                pm.LeaveOpen = true;
                if nargin < 2
                    seed = 1;
                    if nargin < 1
                        oga = false;
                    end
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            
            dim = 1;
            x = repmat(-10:.1:10,dim,1);
            [m,M] = general.Utils.getBoundingBox(x);
            dia = Norm.L2(M-m)/10;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.NumGammas = 1;
            a.Dists = dia;
            a.MaxExpansionSize = 100;
            a.UsefScaling = false;
            a.UseOGA = oga;
            a.MaxRelErr = 1e-2;
            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            % Use same dist aka space
            f.Kernel.setGammaForDistance(a.Dists,a.gameps);
            f.Centers.xi = r.rand(dim,6)*20-10;
            f.Ma = r.rand(dim,6)*5;%-2.5;
            M = f.MBnd;
            
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
            
            x = repmat(-15:.01:15,dim,1);
            fx = f.evaluate(x);
            afx = kexp.evaluate(x);
            err2 = Norm.L2((fx-afx)');
            pm.nextPlot('vkoga_approx');
            plot(x',fx','r',x',afx','b',...
                f.Centers.xi, f.evaluate(f.Centers.xi), 'k.',...
                kexp.Centers.xi,f.evaluate(kexp.Centers.xi),'rx','MarkerSize',10);
            if dim == 1
                legend('f(x)','kexp(x)','f(x) Centers','Centers at kexp(x)');
            end
            title(sprintf('Stopping condition: Max L^2-error on training data: %e\nTotal L^2-error: %e',...
                a.MaxRelErr, err2));
            
            pm.nextPlot('vkoga_bounds');
            n = a.expsizes(1);
            %bound = M^2 ./ (1:n);
            bound = dim*M^2 ./ (1+(1:n)/dim);
            semilogy(1:n,f.NativeNorm^2 - a.HerrDecay(1:n),'r',1:n, bound,'b');%,1:n, bound2,'b--');
            legend('H-approximation error at step m','Theoretical upper bound');
            pm.done;
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
            %[m,M] = general.Utils.getBoundingBox(x);
%             dia = Norm.L2(M-m);
%             dia = sqrt(dim)*xrange;
            dia = dim*xrange;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            a.Dists = sqrt(dia);
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,200);
            a.UsefScaling = false;
            a.MaxRelErr = 1e-3;
            a.PhiNormMin = sqrt(eps);
            
            vx = unique(r.rand(vxsize,dim)*xrange+xoff,'rows')';
            a.vxtmuargs{1} = vx;
            
            ao = a.clone;
            ao.UseOGA = true;
            
            % Setup test functions            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(sqrt(dia),a.gameps)
            
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
        
        function plotStatistics(a, ao, f, pm)
            n1 = a.expsizes;
            n2 = ao.expsizes;
            rows = 1;
            if ~isempty(a.vxtmuargs)
                rows = 2;
            end
            if nargin < 4
                pm = tools.PlotManager(false, rows, 3);
            end
            
            pm.nextPlot('abs_l2_err_train');
            semilogy(1:n1,a.err(1:n1),'r',1:n2,ao.err(1:n2),'b');
            title('Absolute L^2 errors on training set');
            xlabel('expansion size N'); ylabel('absolute error');
            legend('VKOGA','WSOGA2');

            pm.nextPlot('rel_l2_err_train');
            semilogy(1:n1,a.relerr(1:n1),'r',1:n2,ao.relerr(1:n2),'b');
            %title(sprintf('Max Relative errors on training set\nStopping condition: relerr < %e',a.MaxRelErr));
            title('Relative L^2 errors on training set');
            legend('VKOGA','WSOGA2');
            xlabel('expansion size N'); ylabel('relative error');

            M = f.MBnd;
            pm.nextPlot('H_err_decay_bounds');
            fnsq = f.NativeNorm^2;
            semilogy(0:n1,fnsq-[0 a.HerrDecay(1:n1)],'r',0:n2,fnsq-[0 ao.HerrDecay(1:n2)],'b');
            hold on;
            dim = size(f.Ma,1);
            bound_sq = (dim*M^2) ./ (1 + (0:(max(n1,n2))) / dim);
            semilogy(0:max(n1,n2),bound_sq(1:max(n1,n2)+1),'b--');
            bound_cm = M^2 * [dim a.VKOGABound];
            semilogy(0:n1,bound_cm(1:n1+1),'r--');
            %title(sprintf('H-norm decay plus upper convergence bound with M=%3.2f', M));
            title(sprintf('H-norm decay plus upper convergence bound with M=%3.2f', M));
            legend('VKOGA','WSOGA2','Upper bnd','VKOGA a-post. bnd');
            hold off;
            
            if ~isempty(a.vxtmuargs)
                %vsiz = size(a.vfx,2);
                pm.nextPlot('abs_l2_err_val');
                semilogy(1:n1,a.verr(1:n1),'r',1:n2,ao.verr(1:n2),'b');
                title('Absolute L^2 errors on validation set');
                xlabel('expansion size N'); ylabel('absolute errors');
                legend('VKOGA','WSOGA2');

                pm.nextPlot('rel_l2_err_val');
                semilogy(1:n1,a.vrelerr(1:n1),'r',1:n2,ao.vrelerr(1:n2),'b');      
                title('Relative L^2 errors on validation set');
                legend('VKOGA','WSOGA2');
                xlabel('expansion size N'); ylabel('relative errors');
            end
            pm.done;
        end
        
        function plotVKOGARes(res, pm)
            if nargin < 2
                pm = tools.PlotManager(false, 2, 3);
                pm.FilePrefix = 'plotVKOGARes';
            end
            runs = length(res.terr);
            runx = 1:runs;
            pm.nextPlot('errcomp');
            plot(runx, [res.terr; res.terro; res.verr; res.verro]);
            VKOGA_b = res.verr <= res.verro; 
            hold on;
            legend('VKOGA','WSOGA2','VKOGA (val)','WSOGA2 (val)');
            title(sprintf('Maximum L^2-errors on training/validation set\nVKOGA <= WSOGA2 on validation set: %2.2f%%',...
                100*sum(VKOGA_b)/runs));
            xlabel('test run'); ylabel('absolute L2 errors');
            
            pm.nextPlot('relerrcomp');
            plot(runx, [res.vrelerr; res.vrelerro]);
            VKOGA_b = res.vrelerr <= res.vrelerro;
            hold on;
            legend('VKOGA (val)','WSOGA2 (val)');
            title(sprintf('Maximum relative L^2-errors on validation set\nVKOGA <= WSOGA2: %2.2f%%',...
                100*sum(VKOGA_b)/runs));
            xlabel('test run'); ylabel('relative L2 errors');
            
            pm.nextPlot('expsizes');
            plot(runx, [res.cnum; res.cnumo]);
            VKOGA_b = res.cnum <= res.cnumo;
            hold on;
            plot(runx(VKOGA_b),res.cnum(VKOGA_b),'g.','MarkerSize',20);
            plot(runx(~VKOGA_b),res.cnumo(~VKOGA_b),'r.','MarkerSize',20);
            title(sprintf('Expansion sizes at MaxRelErr=%e\nVKOGA better than WSOGA2: %2.2f%%',...
                res.a.MaxRelErr,100*sum(VKOGA_b)/runs));
            legend('VKOGA','WSOGA2');
            xlabel('test run'); ylabel('expansion sizes');
            
            pm.nextPlot('Herr_and_bounds');
            semilogy(runx, res.ogabnd, 'b--', runx, res.vkogabnd, 'r--');
            hold on;
            semilogy(runx, res.Herro, 'b', runx, res.Herr, 'r');
            VKOGA_b = res.Herr <= res.Herro;
            hold on;
            semilogy(runx(VKOGA_b),res.Herr(VKOGA_b),'g.','MarkerSize',20);
            semilogy(runx(~VKOGA_b),res.Herro(~VKOGA_b),'r.','MarkerSize',20);
            title(sprintf('H-errors and bounds at MaxRelErr=%e\nVKOGA better than WSOGA2: %2.2f%%',...
                res.a.MaxRelErr,100*sum(VKOGA_b)/runs));
            legend('VOGA bound (a-pri)','VKOGA bound (a-post)','H-err WSOGA2','H-err VKOGA');
            xlabel('test run'); ylabel('H-errors');
            
            pm.nextPlot('errbnd_ratios');
            hlp = res.ogabnd ./ res.vkogabnd;
            plot(runx, hlp, 'r', runx, repmat(mean(hlp),1,runs), 'r--');
            hold on;
            hlp = res.Herro ./ res.Herr;
            plot(runx, hlp, 'b', runx, repmat(mean(hlp),1,runs), 'b--');
            legend('WSOGA2 to VKOGA bound ratio','average','WSOGA2 to VKOGA Herr ratio','average');
            title('WSOGA2 to VKOGA ratios');
            
            if isfield(res,'gammas')
                pm.nextPlot('gammas');
                plot(1:length(res.gammas),res.gammas);
                legend('gamma values');
                xlabel('run');
                ylabel('gamma value');
            end
        end
        
        function [kexp, kexpo, a, ao, atd, pm] = test_VKOGA_TestFuns(name, atdsize, seed)
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
            %                 a.vfx = fun(vx);
            %                 ao.vfx = a.vfx;
            
            % Original greedy version
            ao.computeApproximation(kexp, atd);
            kexpo = kexp.clone;
            % Normalized version
            a.computeApproximation(kexp, atd);
            
            save(sprintf('vkoga_testfun_%s',name),'a','ao','kexp','kexpo','atd');
            
            pm = tools.PlotManager(false,2,2);
            pm.FilePrefix = sprintf('testfun_%s',name);
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
            
            pm.done;
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
            a.NumGammas = 20;
            a.MaxExpansionSize = 200;
            a.UsefScaling = false;
            a.MaxRelErr = 1e-3;
            a.PhiNormMin = 1e-3;%sqrt(eps);
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