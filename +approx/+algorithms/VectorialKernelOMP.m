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
        PhiNormMin = eps;
        
        % The original kernel expansion used to generate the approximation
        % training data (DEBUG)
        f;
        
        VKOGABound;
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
            copy.f = this.f;
            copy.VKOGABound = this.VKOGABound;
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
            this.VKOGABound = this.err;
            
            fxinorm = Norm.L2(atd.fxi);
            % Validation set initializations
            if ~isempty(this.vxtmuargs)
                this.verr = this.err;
                this.vrelerr = this.err;
                vfxnorm = Norm.L2(this.vfx);
            end
            
            if ~this.UseOGA && ~isempty(this.f)
                hlp1 = this.f.NativeNorm^2;
                hlp2 = size(this.f.Ma,1)*this.f.M^2;
                fprintf('Initial condition ||f||^2 (%e) <= qM^2 (%e)\n', hlp1, hlp2);
                if hlp1 > hlp2
                    error('Initial condition violated.');
                end
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
                    fDotPhiSqAll_comp = (atd.fxi(:,free) - F*A).^2;
                    
                     % Sum
                    if this.UseOGA
                        fDotPhiSqAll = sum(fDotPhiSqAll_comp,1);
                    else
                        fDotPhiSqAll = sum(fDotPhiSqAll_comp,1) ./ phinormsq;
                    end
                    
                    if ~this.UseOGA && siz > 0 && ~isempty(this.f) && false
                        M = this.f.M;
                        % M-estimation verification stuff
                        f_fm = this.f - kexp;
%                         cn = f_fm.ComponentNorms.^2;
                        cn = sum(this.f.Ma .* (this.f.evaluate(this.f.Centers.xi) - kexp.evaluate(this.f.Centers.xi)),2);
%                         s3 = sum(abs(this.f.Ma) .* abs((this.f.evaluate(this.f.Centers.xi) - kexp.evaluate(this.f.Centers.xi))),2);
%                         s3b = sum(abs(this.f.Ma),2) .* max(abs((this.f.evaluate(this.f.Centers.xi) - kexp.evaluate(this.f.Centers.xi))),[],2);
%                         s3c = sum(abs(this.f.Ma),2) .* max(abs((this.f.evaluate(atd.xi) - kexp.evaluate(atd.xi))),[],2);
%                         A = i.interpolate(kexp.getKernelVector(atd.xi));
%                         F = atd.fxi(:,used);
%                         fdotphi = abs(atd.fxi - F*A);
%                         s4b = sum(abs(this.f.Ma),2) .* max(fdotphi,[],2);
%                         s4 = sum(abs(this.f.Ma),2) .* max(sqrt(fDotPhiSqAll_comp),[],2);
                        ubnd = M * max(sqrt(fDotPhiSqAll_comp),[],2);
                        pt = PrintTable;
                        pt.addRow('||f_j-f_j^m||^2', 'M*max <f_j,phi~>', 'diff');
                        pt.addRow(cn, ubnd, ubnd-cn);
%                         pt.display;
                        if any(ubnd < cn)
                            warning('fck:id','Upper bound estimation for gain failed.');
                        end
%                         
%                         % sum/max exchange verification stuff
                        dim = size(this.f.Ma,1);
                        maxsum = max(sum(fDotPhiSqAll_comp,1));
                        summax = sum(max(fDotPhiSqAll_comp,[],2));
                        pt = PrintTable;
                        pt.addRow('Dim', 'Gain','MaxSum','SumMax','sum(f_j-f_j^m)^4/(qM^2)','sum(cn)^2/(d*M)^2','f_fm.NativeNorm^4/(d*M)^2');
                        hlp = f_fm.NativeNorm^4/(dim*M)^2;
                        pt.addRow(dim, v, maxsum, summax/dim, 1/(dim*M^2) * sum(cn.^2), ...
                            sum(cn)^2/(dim*M)^2, hlp, {'%2.3e'});
                        pt.ColSep = ' >= ';
%                         pt.display;
                        if v < hlp
                            error('Gain term smaller than lower bound');
                        end
                        
                        fd = this.f - kexp;
                        fdprev = this.f - kexp_prev;
                        pt = PrintTable;
                        pt.addRow('||f-f^m||^2','||f-f^{m-1}||^2 - gain','||f-f^{m-1}||^2(1-||f-f^{m-1}||^2/q^2M^2)', '||f-f^{m-1}||^2','gain');
                        fdprevn2 = fdprev.NativeNorm^2;
                        hlp = fdprevn2*(1-fdprevn2/(dim*M)^2);
                        pt.addRow(fd.NativeNorm^2, fdprevn2 - v,...
                            hlp, fdprevn2, v, {'%2.3e'});
                        pt.ColSep = ' <= ';
                        pt.display;
                        if fd.NativeNorm^2 > hlp || fdprevn2 - v > hlp
                            error('Wrong gain calculations');
                        end
                        
                        pt = PrintTable;
                        pt.addRow('||f-f^m||^2','||f||^2 - sum_k^m sum_j^d<f_j,phi_x^{m-1}>^2');
                        hlp1 = fd.NativeNorm^2;
                        hlp2 = this.f.NativeNorm^2 - this.HerrDecay(gidx, siz);
                        pt.addRow(hlp1, hlp2);
                        pt.display;
                        if abs(hlp1-hlp2)/abs(hlp1) > 5e-4
                            warning('freak:out','Inconsistent ||f-f^m|| vs. ||f|| - accumulated gain');
                        end
                    end
                    
                    [v, maxidx] = max(fDotPhiSqAll);
                    
                    used(end+1) = free(maxidx);%#ok
                    kexp_prev = kexp.clone;
                    this.extendExpansion(kexp, atd, used(end));
                    siz = siz+1;
                    
                    if this.UseOGA
                        this.Gain(gidx,siz) = v / phinormsq(maxidx);
                    else
                        this.Gain(gidx,siz) = v;
                    end
                    off = 0; off2 = 0;
                    if siz > 1
                        off = this.HerrDecay(siz-1);
                        off2 = this.VKOGABound(siz-1);
                    end
                    this.HerrDecay(gidx, siz) = off + this.Gain(gidx,siz);
                    % Only build improved bound if VKOGA is used
                    if ~this.UseOGA
                        this.VKOGABound(gidx, siz) = off2 + 1/max(phinormsq);
                    end
                    
                    %% Debug - change later to better stopping criteria not using the kernel expansion
                    i.init(kexp.getKernelMatrix);
                    kexp.Ma = i.interpolate(atd.fxi(:,used))';
                    
                    if siz > 1 %&& KerMor.App.Verbose > 1
                        if KerMor.App.Verbose > 3
                            this.doplots(Kbig, A, kexp, afxi, phinormsq, atd, ...
                                free, fDotPhiSqAll, maxidx, v, gidx, siz);
                        end
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
            % Clean up: remove improved bound if algorithm does use classic
            % OGA
            if ~this.UseOGA
                dim = size(kexp.Ma,1);
                this.VKOGABound = dim ./ (1 + this.VKOGABound/dim);
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