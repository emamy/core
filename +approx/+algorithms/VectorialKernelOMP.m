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
    end
    
    methods
        function this = VectorialKernelOMP
            this = this@approx.algorithms.BaseAdaptiveCWKA;
        end

        function copy = clone(this)
            % Clones the instance.
            copy = approx.algorithms.VectorialKernelOMP;
            copy = clone@approx.algorithms.BaseAdaptiveCWKA(this, copy);
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedAdaptiveApproximation(this, kexp, atd)
            
            % Sum
            [~, initialidx] = max(sum(abs(atd.fxi),1));
            % Supremum
            %[v, idx] = max(max(abs(atd.fxi),[],1));
            
            i = general.interpolation.KernelInterpol;
            total = 1:size(atd.xi,2);
            
            if isempty(this.Dists)
                d = this.getDists(atd);
            else
                d = this.Dists;
            end
            this.err = zeros(this.NumGammas,this.MaxExpansionSize);
            this.relerr = this.err;
            this.Gain = this.err;
            % Make one bigger as it starts with zero gain
            this.HerrDecay = [zeros(this.NumGammas,1) this.Gain];
            fxinorm = sqrt(sum(atd.fxi.^2));
            
            minHerr = Inf;
            for gidx = 1:size(d,2)
                
                this.setDistKernelConfig(kexp, d(:,gidx));
                
                % Clear kexp
                kexp.Centers.xi = [];
                if this.pte
                    kexp.Centers.ti = [];
                    kexp.Centers.mui = [];
                end
                kexp.Ma = [];
                % Add first vector.
                this.extendExpansion(kexp, atd, initialidx);
                used = initialidx;
                
                %% Debug
                K = kexp.getKernelMatrix;
                i.init(K);
                kexp.Ma = i.interpolate(atd.fxi(:,used))';
                afxi = kexp.evaluate(atd.xi);
                e = sqrt(sum((afxi-atd.fxi).^2));
                this.err(gidx,1) = max(e);
                this.relerr(gidx,1) = max(e./fxinorm);
                projpart_old = -sum(sum(kexp.Ma .* atd.fxi(:,used),2));
                %projpart_old = projpart_old + kexp.Ma'*(K*kexp.Ma);
                
                siz = 1;
                while siz < this.MaxExpansionSize && this.relerr(gidx,siz) > this.MaxRelErr
                    % get indices of still free centers
                    free = setdiff(total,used);
                    
                    % Compute LaGrange basis coefficients
                    K = kexp.getKernelMatrix;
                    %                 if siz > 1
                    %                     mesh(K);
                    %                 end
                    
                    %                 [P, k2] = i.getPreconditioner(kexp.Kernel, kexp.Centers.xi);
                    %                 fprintf('cond(K)=%e, cond(P*K)=%e, k2=%d\n',cond(K),cond(P*K),k2);
                    %                 K = P*K; I = P*I;
                    
                    args = {};
                    if this.pte
                        args = {atd.ti(free), atd.mui(:,free)};
                    end
                    Kbig = kexp.getKernelVector(atd.xi(:,free),args{:})';
                    
                    i.init(data.MemoryKernelMatrix(K));
                    
                    % Compute kernel fcn projection to H(m-1)
                    A = i.interpolate(Kbig');
                    F = atd.fxi(:,used);
                    
                    % Squared versions
                    fDotPhiSqAll = (atd.fxi(:,free) - F*A).^2;
                    phinormsq = 1 - sum(A.*Kbig,1);
                    % normal versions
                    %                     fDotPhiSqAll = abs(atd.fxi(:,free) - F*A);
                    %                     phinormsq = sqrt(1 - sum(A.*Kbig,1));
                    % same result, but which is numerically cleverer?
                    
                    % Cap too small norms!
                    %phinormsq = max(phinormsq,1e-2);
                    
                     % Sum
                    if this.UseOGA
                        fDotPhiSqAll = sum(fDotPhiSqAll,1);
                    else
                        fDotPhiSqAll = sum(fDotPhiSqAll,1) ./ phinormsq;
                    end
                    
                    [v, maxidx] = max(fDotPhiSqAll);
                    
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
                    
                    if KerMor.App.Verbose > 3
                        this.doplots(Kbig, A, kexp, afxi, phinormsq, atd, ...
                            free, fDotPhiSqAll, maxidx, v, gidx, siz);
                    end
                    
                    used(end+1) = free(maxidx);%#ok
                    this.extendExpansion(kexp, atd, used(end));
                    
                    siz = siz+1;
                    
                    %% Debug - change later to better stopping criteria not using the kernel expansion
                    K = kexp.getKernelMatrix;
                    i.init(K);
                    kexp.Ma = i.interpolate(atd.fxi(:,used))';
                    %                     kexp.Ma = (kexp.getKernelMatrix \ atd.fxi(:,used)')';
                    afxi = kexp.evaluate(atd.xi);
                    e = sqrt(sum((afxi-atd.fxi).^2));
                    this.err(gidx,siz) = max(e);
                    this.relerr(gidx,siz) = max(e./fxinorm);
                    
                    projpart = -sum(sum(kexp.Ma .* atd.fxi(:,used),2));
                    %projpart = projpart + kexp.Ma*(K*kexp.Ma');
                    
                    gain_direct = abs(projpart - projpart_old);
                    gain_impl = this.Gain(gidx,siz-1);
                    fprintf('Direct gain: %e, implicit gain: %e, difference: %e, ||f^m||_H^2: %e\n',...
                        gain_direct,gain_impl,gain_direct-gain_impl, kexp.NativeNorm^2);
                    
                    projpart_old = projpart;
                end
                siz = siz-1;
                
                if this.err(gidx,siz) < minHerr
                    bestused = used;
                    bestMa = kexp.Ma;
                    bestDist = d(:,gidx);
                    minHerr = this.err(gidx,siz);
                    bestgidx = gidx;
                end
            end
            
            % Remove first zero column
            this.HerrDecay(:,1) = [];
            %this.HerrDecay = repmat(max(this.HerrDecay,[],2),1,this.MaxExpansionSize) - this.HerrDecay;
            
            g = this.setDistKernelConfig(kexp, bestDist);
            og = 'off';
            if this.UseOGA
                og = 'on';
            end
            if KerMor.App.Verbose > 1
                fprintf('VKOMP best gamma index:%d (value:%e) for VKOMP with OGA=%s, exp-size:%d\n',...
                    bestgidx,g,og,size(bestMa,2));
            end
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
            cafxi = kexp.evaluate(cxi);
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
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            %a.CoeffComp.UseLU = true;
            a.NumGammas = 1;
            a.Dists = 2;
            a.MaxExpansionSize = 10;
            %a.ValidationPercent = .0001;
            a.UsefScaling = false;
            a.UseOGA = oga;
            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            % Use same dist aka space
            f.Kernel.setGammaForDistance(a.Dists,a.gameps);
            f.Centers.xi = r.rand(1,6)*20-10;
            f.Ma = r.rand(1,6)*5;%-2.5;
            M = sum(abs(f.Ma));
            
            x = -10:.1:10;
            fx = f.evaluate(x);
            atd = data.ApproxTrainData(x,[],[]);
            atd.fxi = fx;
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            a.computeApproximation(kexp, atd);
            
            x = -20:.01:20;
            fx = f.evaluate(x);
            afx = kexp.evaluate(x);
            err2 = sqrt(sum((fx-afx).^2));
            figure(3);
            subplot(1,2,1);
            plot(x,fx,'r',x,afx,'b',...
                f.Centers.xi, f.evaluate(f.Centers.xi), 'k.',...
                kexp.Centers.xi,f.evaluate(kexp.Centers.xi),'rx','MarkerSize',10);
            legend('f(x)','kexp(x)','f(x) Centers','Centers at kexp(x)');
            title(sprintf('Total L^2-error: %e',err2));
            axis tight;
            subplot(1,2,2);
            n = length(a.HerrDecay);
            plot(1:n,f.NativeNorm^2 - a.HerrDecay,'r',1:n, M^2 ./ (1:n),'b');
            legend('Approximation error at step m in H','Theoretical upper bound M/sqrt(m)');
            axis tight;
        end
        
        function [res, kexp, kexp_OGA, a, f, atd] = test_VKOMP_Versions(runs, seed)
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
            centers = 10;
            dim = 1;
            atdsize = 100;

            x = unique(r.rand(atdsize,dim)*10-5,'rows')';
            [m,M] = general.Utils.getBoundingBox(x);
            dia = sqrt(sum((M-m).^2));
            
            a = approx.algorithms.VectorialKernelOMP;
            a.CoeffComp = general.interpolation.KernelInterpol;
            %a.CoeffComp.UseLU = true;
            a.Dists = sqrt(dia);
            a.NumGammas = 1;
            a.MaxExpansionSize = 500;
            a.UsefScaling = false;
            a.MaxRelErr = 1e-4;
            
            % Setup test functions            
            f = kernels.KernelExpansion;
            f.Kernel = kernels.GaussKernel;
            f.Kernel.setGammaForDistance(sqrt(dia),a.gameps);
            atd = data.ApproxTrainData(x,[],[]);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel;
            
            vx = unique(r.rand(2000,dim)*20-10,'rows')';
            res = zeros(8,runs);
            pi = tools.ProcessIndicator('Running tests', runs, false);
            for run=1:runs
                f.Centers.xi = r.rand(dim,centers)*20-10;
                f.Ma = r.rand(dim,centers)*15;%-2.5;
                atd.fxi = f.evaluate(x);
                
                % Evaluate error on validation set
                fx = f.evaluate(vx);
                fxinorm = sqrt(sum(fx.^2));
                
                % Original greedy version
                a.UseOGA = true;
                a.computeApproximation(kexp, atd);
                res(2,run) = a.err(size(kexp.Ma,2));
                afx = kexp.evaluate(vx);
                e = sqrt(sum((fx-afx).^2));
                res(4,run) = max(e);
                res(8,run) = max(e./fxinorm);
                res(6,run) = size(kexp.Ma,2);
                
                kexp_OGA = kexp.clone;
                oerr = a.err;
                orelerr = a.err;
                odec = a.HerrDecay;
                
                % Normalized version
                a.UseOGA = false;
                a.computeApproximation(kexp, atd);
                afx = kexp.evaluate(vx);
                res(1,run) = a.err(size(kexp.Ma,2));
                e = sqrt(sum((fx-afx).^2));
                res(3,run) = max(e);
                res(7,run) = max(e./fxinorm);
                res(5,run) = size(kexp.Ma,2);
                
                pi.step(run);
                
                if runs == 1
                    n1 = length(a.err);
                    n2 = length(oerr);
                    figure(12);
                    subplot(1,3,1);
                    semilogy(1:n1,a.err,'r',1:n2,oerr,'b');
                    title('Max Absolute errors on training set');
                    legend('KOMP','OGA');
                    
                    subplot(1,3,2);
                    semilogy(1:n1,a.relerr,'r',1:n2,orelerr,'b');
                    title('Max Relative errors on training set');
                    legend('KOMP','OGA');
                    
                    M = max(sum(abs(f.Ma),2));
                    subplot(1,3,3);
                    
                    plot(1:n1,M^2-a.HerrDecay,'r',1:n2,M^2-odec,'b');
                    hold on;
                    bound_sq = (dim*M^2) ./ (1+ (1:(max(n1,n2))) / dim);
                    plot(1:max(n1,n2),bound_sq(1:n1),'r--');
                    title(sprintf('Decay plus upper convergence bound with M=%3.2f',M));
                    legend('KOMP','OGA','Upper bound');
                end
            end
            pi.stop;
            
            figure;
            subplot(1,3,1);
            plot(1:runs, res(1:4,:));
            komp_b = res(3,:) <= res(4,:); 
            hold on;
            legend('KOMP','OGA','KOMP (val)','OGA (val)');
            title(sprintf('Maximum L^2-errors on training/validation set\nKOMP <= OGA on validation set: %2.2f%%',...
                100*sum(komp_b)/runs));
            xlabel('test run'); ylabel('absolute L2 errors');
            axis tight;
            
            subplot(1,3,2);
            plot(1:runs, res(7:8,:));
            komp_b = res(7,:) <= res(8,:);
            hold on;
            legend('KOMP (val)','OGA (val)');
            title(sprintf('Maximum relative L^2-errors on validation set\nKOMP <= OGA: %2.2f%%',...
                100*sum(komp_b)/runs));
            xlabel('test run'); ylabel('relative L2 errors');
            axis tight;
            
            runx = 1:runs;
            subplot(1,3,3);
            plot(runx,res(5:6,:));
            komp_b = res(5,:)-res(6,:) <= 0;
            hold on;
            plot(runx(komp_b),res(5,komp_b),'g.','MarkerSize',20);
            plot(runx(~komp_b),res(6,~komp_b),'r.','MarkerSize',20);
            title(sprintf('Expansion sizes at MaxRelErr=%e\nKOMP better than OGA: %2.2f%%',...
                a.MaxRelErr,100*sum(komp_b)/runs));
            legend('KOMP','OGA');
            xlabel('test run'); ylabel('expansion sizes');
            axis tight;
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