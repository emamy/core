classdef VKOMP
% VKOMP: 
%
% Tests for the VKOGA algorithm
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
        function [kexp, atd, a, f] = test_VKOMP(oga, seed)
            if nargin < 2
                seed = 1;
                if nargin < 1
                    oga = false;
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            
            dim = 2;
            x = repmat(-10:.1:10,dim,1);
            [m,M] = general.Utils.getBoundingBox(x);
            dia = Norm.L2(M-m)/10;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.NumGammas = 1;
            a.Dists = dia;
            a.MaxExpansionSize = 300;
            a.UsefScaling = false;
            a.UseOGA = oga;
            a.MaxRelErr = 1e-3;
            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            % Use same dist aka space
            f.Kernel.setGammaForDistance(a.Dists,a.gameps);
            f.Centers.xi = r.rand(dim,6)*20-10;
            f.Ma = r.rand(dim,6)*5;%-2.5;
            M = max(Norm.L1(f.Ma'));
            
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
            
            x = repmat(-20:.01:20,dim,1);
            fx = f.evaluate(x);
            afx = kexp.evaluate(x);
            err2 = sum(Norm.L2(fx-afx));
            figure(3);
            subplot(1,2,1);
            plot(x',fx','r',x',afx','b',...
                f.Centers.xi, f.evaluate(f.Centers.xi), 'k.',...
                kexp.Centers.xi,f.evaluate(kexp.Centers.xi),'rx','MarkerSize',10);
            if dim == 1
                legend('f(x)','kexp(x)','f(x) Centers','Centers at kexp(x)');
            end
            title(sprintf('Total L^2-error: %e',err2));
            axis tight;
            subplot(1,2,2);
            n = a.expsizes(1);
            bound = M^2 ./ (1:n);
            bound2 = dim*M^2 ./ (1+(1:n)/dim);
            plot(1:n,f.NativeNorm^2 - a.HerrDecay(1:n),'r',1:n, bound,'b',1:n, bound2,'b--');
            legend('Approximation error at step m in H','Theoretical upper bound M/sqrt(m)');
            axis tight;
        end
        
        function [res, kexp, kexp_OGA, a, ao, f, atd] = test_VKOMP_Versions(runs, seed)
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
            dim = 100;
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
            a.PhiNormMin = eps;
            
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
                s = size(kexp.Ma,2);
                res.terro(:,run) = ao.err(:,s);
                res.verro(:,run) = ao.verr(:,s);
                res.vrelerro(:,run) = ao.vrelerr(:,s);
                res.cnumo(run) = size(kexp.Ma,2);
                
                kexp_OGA = kexp.clone;
                
                % Normalized version
                a.computeApproximation(kexp, atd);
                s = size(kexp.Ma,2);
                res.terr(:,run) = a.err(:,s);
                res.verr(:,run) = a.verr(:,s);
                res.vrelerr(:,run) = a.vrelerr(:,s);
                res.cnum(run) = size(kexp.Ma,2);
                
                pi.step(run);
                
                if runs == 1 || true
                    testing.VKOMP.plotStatistics(a,ao,f);
                    if runs > 1
                        pause;
                    end
                end
            end
            pi.stop;
            
            if nargout == 0
                testing.VKOMP.plotKompRes(res);
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
        
        function plotStatistics(a, ao, f)
            persistent h;
            n1 = a.expsizes;
            n2 = ao.expsizes;
            rows = 1;
            if ~isempty(a.vxtmuargs)
                rows = 2;
            end
            
            if isempty(h) || ~ishandle(h)
                h = figure;
            else
                figure(h);
            end
            subplot(rows,3,1);
            semilogy(1:n1,a.err(1:n1),'r',1:n2,ao.err(1:n2),'b');
            title('Max Absolute errors on training set');
            xlabel('Expansion size'); ylabel('Absolute L2 error');
            legend('KOMP','OGA'); axis tight;

            subplot(rows,3,2);
            semilogy(1:n1,a.relerr(1:n1),'r',1:n2,ao.relerr(1:n2),'b');
            title(sprintf('Max Relative errors on training set\nStopping condition: relerr < %e',a.MaxRelErr));
            legend('KOMP','OGA'); axis tight;
            xlabel('Expansion size'); ylabel('Relative L2 error');

            M = f.M;
            subplot(rows,3,3);
            fnsq = f.NativeNorm^2;
            semilogy(0:n1,fnsq-[0 a.HerrDecay(1:n1)],'r',0:n2,fnsq-[0 ao.HerrDecay(1:n2)],'b');
            hold on;
            dim = size(f.Ma,1);
            bound_sq = (dim*M^2) ./ (1 + (0:(max(n1,n2))) / dim);
            semilogy(0:max(n1,n2),bound_sq(1:max(n1,n2)+1),'b--');
            psum = zeros(1,n1);
            psum(1) = 1/a.phinormmax(1);
            for i=2:n1
                psum(i) = psum(i-1) + 1/a.phinormmax(i);
            end
            bound_cm = (dim*M^2) ./ (1 + [0 psum]/dim);
            semilogy(0:n1,bound_cm,'r--');
            title(sprintf('Decay plus upper convergence bound with M=%3.2f', M));
            legend('KOMP','OGA','Upper bound','VKOGA a-post bnd');
            axis tight;
            hold off;
            
            if ~isempty(a.vxtmuargs)
                vsiz = size(a.vfx,2);
                subplot(2,3,4);
                semilogy(1:n1,a.verr(1:n1),'r',1:n2,ao.verr(1:n2),'b');
                title(sprintf('Max Absolute errors on validation set (size %d)',vsiz));
                axis tight;
                xlabel('Expansion size'); ylabel('Absolute L2 error');
                legend('KOMP','OGA');

                subplot(2,3,5);
                semilogy(1:n1,a.vrelerr(1:n1),'r',1:n2,ao.vrelerr(1:n2),'b');      
                title(sprintf('Max Relative errors on validation set (size %d)',vsiz));
                legend('KOMP','OGA');
                xlabel('Expansion size'); ylabel('Relative L2 error');
                axis tight;
            end
        end
        
        function plotKompRes(res)
            figure;
            runs = length(res.terr);
            runx = 1:runs;
            subplot(1,3,1);
            plot(runx, [res.terr; res.terro; res.verr; res.verro]);
            komp_b = res.verr <= res.verro; 
            hold on;
            legend('KOMP','OGA','KOMP (val)','OGA (val)');
            title(sprintf('Maximum L^2-errors on training/validation set\nKOMP <= OGA on validation set: %2.2f%%',...
                100*sum(komp_b)/runs));
            xlabel('test run'); ylabel('absolute L2 errors');
            axis tight;
            
            subplot(1,3,2);
            plot(runx, [res.vrelerr; res.vrelerro]);
            komp_b = res.vrelerr <= res.vrelerro;
            hold on;
            legend('KOMP (val)','OGA (val)');
            title(sprintf('Maximum relative L^2-errors on validation set\nKOMP <= OGA: %2.2f%%',...
                100*sum(komp_b)/runs));
            xlabel('test run'); ylabel('relative L2 errors');
            axis tight;
            
            subplot(1,3,3);
            plot(runx, [res.cnum; res.cnumo]);
            komp_b = res.cnum <= res.cnumo;
            hold on;
            plot(runx(komp_b),res.cnum(komp_b),'g.','MarkerSize',20);
            plot(runx(~komp_b),res.cnumo(~komp_b),'r.','MarkerSize',20);
            title(sprintf('Expansion sizes at MaxRelErr=%e\nKOMP better than OGA: %2.2f%%',...
                res.a.MaxRelErr,100*sum(komp_b)/runs));
            legend('KOMP','OGA');
            xlabel('test run'); ylabel('expansion sizes');
            axis tight;
        end
    end
    
end