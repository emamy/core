classdef VectorialKernelOMP < approx.algorithms.BaseAdaptiveCWKA
% VectorialKernelOMP: 
%
%
%
% @author Daniel Wirtz @date 2012-02-09
%
% @new{0,6,dw,2012-02-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        err;
        relerr;
        
        Gain;
        UseOGA = false;
        HerrDecay;
        Dists;
        
        vxtmuargs;
        vfx;
        verr;
        vrelerr;
        
        used;
        expsizes;
        
        % Lower bound for orthonormal remainders (kernel matrix conditioning improvement for
        % larger bounds
        PhiNormMin = 0;
    end
    
    methods
        function this = VectorialKernelOMP
            this = this@approx.algorithms.BaseAdaptiveCWKA;
        end

        function copy = clone(this)
            % Clones the instance.
            copy = approx.algorithms.VectorialKernelOMP;
            copy = clone@approx.algorithms.BaseAdaptiveCWKA(this, copy);
            copy.PhiNormMin = this.PhiNormMin;
            copy.vxtmuargs = this.vxtmuargs;
            copy.vfx = this.vfx;
            copy.verr = this.verr;
            copy.vrelerr = this.vrelerr;
            copy.Gain = this.Gain;
            copy.UseOGA = this.UseOGA;
            copy.HerrDecay = this.HerrDecay;
            copy.Dists = this.Dists;
            copy.err = this.err;
            copy.relerr = this.relerr;
            copy.used = this.used;
            copy.expsizes = this.expsizes;
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedAdaptiveApproximation(this, kexp, atd)
            
            i = general.interpolation.KernelInterpol;
            total = 1:size(atd.xi,2);
            
            if isempty(this.Dists)
                d = this.getDists(atd);
            else
                d = this.Dists;
            end
            
            %% Debug/error information inits
            xtmuargs = {atd.xi};
            if ~isempty(atd.ti)
                xtmuargs{2} = atd.ti;
            end
            if ~isempty(atd.mui)
                xtmuargs{3} = atd.mui;
            end
            this.err = zeros(this.NumGammas,this.MaxExpansionSize);
            this.relerr = this.err;
            this.Gain = this.err;
            % Make one bigger as it starts with zero gain
            this.HerrDecay = this.err;
            
            fxinorm = Norm.L2(atd.fxi);
            % Validation set initializations
            if ~isempty(this.vxtmuargs)
                this.verr = this.err;
                this.vrelerr = this.err;
                vfxnorm = Norm.L2(this.vfx);
            end
            
            minHerr = Inf;
            %% Run loop for all desired distances
            for gidx = 1:size(d,2)
                
                % Clear kexp
                kexp.clear;
                % Set current hyperconfiguration
                this.setDistKernelConfig(kexp, d(:,gidx));
                
                projpart_old = 0;
                used = []; %#ok<*PROP>
                siz = 0;
                %% Main extension loop
                while true
                    % get indices of still free centers
                    free = total;
                    free(used) = [];
                    
                    args = {atd.xi(:,free)};
                    if ~isempty(atd.ti)
                        args{2} = atd.ti(free);
                    end
                    if ~isempty(atd.mui)
                        args{3} = atd.mui(:,free);
                    end
                    
                    % We only have remainders as of the second run!
                    if siz > 0
                        Kbig = kexp.getKernelVector(args{:})';

                        % Compute kernel fcn projection to H(m-1)
                        i.init(kexp.getKernelMatrix);
                        A = i.interpolate(Kbig');
                        F = atd.fxi(:,used);
                        
                        % Cap too small norms!
                        phinormsq = max(1 - sum(A.*Kbig,1), this.PhiNormMin^2);
                    else
                        F = 0; A = 0; phinormsq = ones(size(total));
                    end
                    
                    % Squared versions
                    fDotPhiSqAll = (atd.fxi(:,free) - F*A).^2;
                    
                     % Sum
                    if this.UseOGA
                        fDotPhiSqAll = sum(fDotPhiSqAll,1);
                    else
                        fDotPhiSqAll = sum(fDotPhiSqAll,1) ./ phinormsq;
                    end
                    
                    [v, maxidx] = max(fDotPhiSqAll);
                    
                    used(end+1) = free(maxidx);%#ok
                    this.extendExpansion(kexp, atd, used(end));
                    siz = siz+1;
                    
                    if this.UseOGA
                        this.Gain(gidx,siz) = v / phinormsq(maxidx);
                    else
                        this.Gain(gidx,siz) = v;
                    end
                    off = 0;
                    if siz > 1
                        off = this.HerrDecay(siz-1);
                    end
                    this.HerrDecay(gidx, siz) = off + this.Gain(gidx,siz);
                    
                    %% Debug - change later to better stopping criteria not using the kernel expansion
                    i.init(kexp.getKernelMatrix);
                    kexp.Ma = i.interpolate(atd.fxi(:,used))';
                    
                    if KerMor.App.Verbose > 3 && siz > 1
                        this.doplots(Kbig, A, kexp, afxi, phinormsq, atd, ...
                            free, fDotPhiSqAll, maxidx, v, gidx, siz);
                    end
                    
                    % Compute errors
                    compTrainingError;
                    compValidationError;
                    
                    projpart = -sum(sum(kexp.Ma .* atd.fxi(:,used),2));
                    
                    gain_direct = abs(projpart - projpart_old);
                    gain_impl = this.Gain(gidx,siz);
                    if KerMor.App.Verbose > 2
                        fprintf('VKOMP direct gain: %e, implicit gain: %e, difference: %e, ||f^m||_H^2: %e\n',...
                            gain_direct,gain_impl,gain_direct-gain_impl, kexp.NativeNorm^2);
                    end
                    projpart_old = projpart;
                    
                    % Check stopping conditions
                    if siz >= this.MaxExpansionSize || this.relerr(gidx,siz) <= this.MaxRelErr
                        this.HerrDecay(gidx, siz+1:end) = this.HerrDecay(gidx, siz);
                        this.err(gidx, siz+1:end) = this.err(gidx, siz);
                        this.relerr(gidx, siz+1:end) = this.relerr(gidx, siz);
                        break;
                    end
                end
                
                this.used(gidx,1:siz) = used;
                this.expsizes(gidx) = siz;
                
                if this.err(gidx,siz) < minHerr
                    bestused = used;
                    bestMa = kexp.Ma;
                    bestDist = d(:,gidx);
                    minHerr = this.err(gidx,siz);
                    bestgidx = gidx;
                end
            end
            
            % Restore best configuration
            g = this.setDistKernelConfig(kexp, bestDist);
            kexp.Ma = bestMa;
            kexp.Centers.xi = atd.xi(:,bestused);
            if this.pte
                if atd.hasTime
                    kexp.Centers.ti = atd.ti(:,bestused);
                end
                if atd.hasParams
                    kexp.Centers.mui = atd.mui(:,bestused);
                end
            end
            
            if KerMor.App.Verbose > 1
                og = 'off';
                if this.UseOGA
                    og = 'on';
                end
                fprintf('VKOMP best gamma index:%d (value:%e) for VKOMP with OGA=%s, exp-size:%d\n',...
                    bestgidx,g,og,size(bestMa,2));
            end
            
            function compTrainingError
                afxi = kexp.evaluate(xtmuargs{:});
                e = Norm.L2(afxi-atd.fxi);
%                 this.err(gidx,siz) = max(e);
%                 this.relerr(gidx,siz) = max(e./fxinorm);
                this.err(gidx,siz) = Norm.L2(e');
                this.relerr(gidx,siz) = Norm.L2((e./fxinorm)');
            end
            
            function compValidationError
                if ~isempty(this.vxtmuargs)
                    avfx = kexp.evaluate(this.vxtmuargs{:});
                    e = Norm.L2(avfx-this.vfx);
%                     this.verr(gidx,siz) = max(e);
%                     this.vrelerr(gidx,siz) = max(e./vfxnorm);
                    this.verr(gidx,siz) = Norm.L2(e');
                    this.vrelerr(gidx,siz) = Norm.L2((e./vfxnorm)');
                end
            end
        end
    end
    
    methods(Access=private)
        function doplots(this, Kbig, A, kexp, afxi, phinormsq, atd, ...
                free, fDotPhiSqAll, maxidx, v, gidx, siz)
            xplotdim = 1;
            yplotdim = 1;
            figure(11);
            subplot(2,3,1);
            plot(Kbig');
            title('Similarity values for new candidates (Kbig)');
            %                         semilogy(this.err(gidx,:));
            %                         title(sprintf('L^2 error on all %d samples, current: %e',size(atd.xi,2),this.err(gidx,siz)));
            axis tight;
            subplot(2,3,2);
            plot(A');
            title('Coefficients A');
            axis tight;
            subplot(2,3,3);
            semilogy(phinormsq);
            legend('phinormsq');
            title('Norms of orthogonalized new basis candidates');
            axis tight;

            pxi = atd.xi(xplotdim,:);
            pxifree = pxi(free);
            cxi = kexp.Centers.xi;
            args = {};
            if this.pte
                args = {kexp.Centers.ti, kexp.Centers.mui};
            end
            cafxi = kexp.evaluate(cxi,args{:});
            cxi = cxi(xplotdim,:);
            cafxi = cafxi(yplotdim,:);

            subplot(2,3,4);
            plot(pxifree,fDotPhiSqAll,'r',pxi,atd.fxi(yplotdim,:),'b',...
                pxi,afxi(yplotdim,:),'g',pxifree(maxidx),v,'b.',cxi,cafxi,'g.');
            title(sprintf('Extension choice %d',free(maxidx)));
            lbl = '<f,phi>^2';
            if this.UseOGA
                lbl = '<f,phi~>^2';
            end
            legend(lbl,'f(x)','approx');
            axis tight;

            subplot(2,3,5);
            semilogy(this.Gain(gidx,1:siz));
            title(sprintf('Gain from %d to %d basis functions: %e',...
                siz-1,siz,this.Gain(gidx,siz)));
            axis tight;

            subplot(2,3,6);
            semilogy(this.HerrDecay(gidx,1:siz));
            title(sprintf('Accumulated H projection error at %d basis functions',siz));
            axis tight;
        end
    end
    
    methods(Static)
        function [kexp, atd, a, f] = test_VKOMP(oga, seed)
            if nargin < 2
                seed = 1;
                if nargin < 1
                    oga = false;
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            
            dim = 5;
            x = repmat(-10:.1:10,dim,1);
            [m,M] = general.Utils.getBoundingBox(x);
            dia = Norm.L2(M-m);
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.NumGammas = 1;
            a.Dists = dia;
            a.MaxExpansionSize = 300;
            a.UsefScaling = false;
            a.UseOGA = oga;
            
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
            n = length(a.HerrDecay);
            plot(1:n,f.NativeNorm^2 - a.HerrDecay,'r',1:n, M^2 ./ (1:n),'b');
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
            rf = RandStream('mt19937ar','Seed',seed);
            centers = 10;
            xrange = 10; xoff = -5;
            dim = 3;
            atdsize = 5000;
            vxsize = 2000;

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
            %dia = Norm.L2(M-m);
            dia = sqrt(dim)*xrange;
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            a.UseOGA = false;
            %a.CoeffComp.UseLU = true;
            a.Dists = sqrt(dia);
            a.NumGammas = 1;
            a.MaxExpansionSize = min(atdsize,200);
            a.UsefScaling = false;
            a.MaxRelErr = 1e-2;
            a.PhiNormMin = 1e-5;
            
            vx = unique(r.rand(vxsize,dim)*20-10,'rows')';
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
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                f.Centers.xi = rf.rand(dim,centers)*20-10;
                f.Ma = rf.rand(dim,centers)*15;%-2.5;
                atd.fxi = f.evaluate(x);
                
                % Set validation f-evaluation
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
                    approx.algorithms.VectorialKernelOMP.plotStatistics(a,ao,f);
                    if runs > 1
                        pause;
                    end
                end
            end
            pi.stop;
            
            if nargout == 0
                approx.algorithms.VectorialKernelOMP.plotKompRes(res);
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
            nativenorm
        end
        
        function plotStatistics(a, ao, f)
            n1 = a.expsizes;
            n2 = ao.expsizes;
            rows = 1;
            if ~isempty(a.vxtmuargs)
                rows = 2;
            end
            
            figure;
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

            M = max(Norm.L1(f.Ma'));
            subplot(rows,3,3);

            plot(1:n1,M^2-a.HerrDecay(1:n1),'r',1:n2,M^2-ao.HerrDecay(1:n2),'b');
            hold on;
            dim = size(f.Ma,1);
            bound_sq = (dim*M^2) ./ (1 + (0:(max(n1,n2))) / dim);
            plot(1:max(n1,n2),bound_sq(1:max(n1,n2)),'r--');
            plot(1:max(n1,n2),1./(1:max(n1,n2)));
            title(sprintf('Decay plus upper convergence bound with M=%3.2f',M));
            legend('KOMP','OGA','Upper bound');
            axis tight;
            
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
            runx = 1:length(res.err);
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
                a.MaxRelErr,100*sum(komp_b)/runs));
            legend('KOMP','OGA');
            xlabel('test run'); ylabel('expansion sizes');
            axis tight;
        end
    end
end

% Compute LaGrange basis coefficients
% K = kexp.getKernelMatrix;
% %                 if siz > 1
% %                     mesh(K);
% %                 end
% I = speye(siz);
% 
% %                 [P, k2] = i.getPreconditioner(kexp.Kernel, kexp.Centers.xi);
% %                 fprintf('cond(K)=%e, cond(P*K)=%e, k2=%d\n',cond(K),cond(P*K),k2);
% %                 K = P*K; I = P*I;
% 
% i.init(data.MemoryKernelMatrix(K));
% A = i.interpolate(I);
% 
% F = atd.fxi(:,used)';
% 
% args = {};
% if this.pte
%     args = {atd.ti(free), atd.mui(:,free)};
% end
% Kbig = kexp.getKernelVector(atd.xi(:,free),args{:});
% 
% fDotPhiSqAll = (atd.fxi(:,free)' - Kbig*A'*F).^2;
% 
% % Sum
% [~, maxidx] = max(sum(fDotPhiSqAll,2));
% % Supremum
% %[v, maxidx] = max(max(abs(fDotPhiSqAll),[],2));
% 
% %plot(sum(fDotPhiSqAll,2))