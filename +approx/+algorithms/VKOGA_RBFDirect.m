classdef VKOGA_RBFDirect < approx.algorithms.AAdaptiveBase
% VKOGA: Vectorial kernel orthogonal greedy algorithm
%
%
%
% @author Daniel Wirtz @date 2012-02-09
%
% @change{0,7,dw,2012-11-26} Renamed to VKOGA and starting to build in IClassConfig interfaces
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
        
        expsizes;
        
        % Lower bound for orthonormal remainders (kernel matrix conditioning improvement for
        % larger bounds
        PhiNormMin = sqrt(eps);
        
        % The original kernel expansion used to generate the approximation
        % training data (DEBUG)
        f;
        
        VKOGABound;
    end
    
    methods
        function this = VKOGA_RBFDirect
            this = this@approx.algorithms.AAdaptiveBase;
        end

        function copy = clone(this)
            % Clones the instance.
            copy = approx.algorithms.VKOGA;
            copy = clone@approx.algorithms.AAdaptiveBase(this, copy);
            copy.PhiNormMin = this.PhiNormMin;
            copy.Gain = this.Gain;
            copy.UseOGA = this.UseOGA;
            copy.HerrDecay = this.HerrDecay;
            copy.relerr = this.relerr;
            copy.used = this.used;
            copy.expsizes = this.expsizes;
            copy.f = this.f;
            copy.VKOGABound = this.VKOGABound;
        end
    end
    
    methods(Access=protected, Sealed)
        function startAdaptiveExtension(this, kexp, atd)
            % Starts the adaptive extension of the VKOGA algorithm.
            
            % Flag for experimental mode
            exp_mode = 1 == 0;
            
            ec = this.ExpConfig;
            nc = ec.getNumConfigurations;
            if nc == 0
                error('Need at least one expansion configuration.');
            end
            
            this.Gain = zeros(nc,this.MaxExpansionSize);
            this.HerrDecay = this.Gain;
            this.VKOGABound = this.Gain;
            
            %% Debug/error information inits
            this.relerr = this.Gain;
            this.MaxErrors = this.Gain;
            
            % Validation set initializations
            if exp_mode && ~this.UseOGA && ~isempty(this.f) && isa(this.f, 'kernels.KernelExpansion')
                this.err = this.Gain;
                hlp1 = this.f.NativeNorm^2;
                hlp2 = size(this.f.Ma,1)*this.f.MBnd^2;
                fprintf('Initial condition ||f||^2 (%e) <= qM^2 (%e)\n', hlp1, hlp2);
                if hlp1 > hlp2
                    error('Initial condition violated.');
                end
            end
            
            minerr = Inf;
            atdfxinormsq = this.ErrorFun(atd.fxi.toMemoryMatrix).^2;
            atdfxinormsq(atdfxinormsq == 0) = 1;
            i = general.interpolation.KernelInterpol;
            total = 1:size(atd.xi,2);
            
            %% Run loop for all desired distances
            pi = ProcessIndicator('VKOGA approximation for %d gamma values',nc,false,nc);
            for cidx = 1:nc
                
                % Re-init expansion from second run on
                if cidx > 1
                    this.initExpansion(kexp, atd);
                end
                siz = 1;
                lastwarn(''); warncnt = 0;
                
                % Set current hyperconfiguration
                ec.applyConfiguration(cidx, kexp);
                
                if exp_mode
                    projpart_old = 0;
                    kexp_prev = kexp;
                end

                %% Main extension loop
                while siz < this.MaxExpansionSize && warncnt < 2
                    
                    % get indices of still free centers
                    free = total;
                    free(this.Used) = [];
                    
                    args = {atd.xi(:,free)};
                    if ~isempty(atd.ti)
                        args{2} = atd.ti(free);
                    end
                    if ~isempty(atd.mui)
                        args{3} = atd.mui(:,free);
                    end
                    
                    Kbig = kexp.getKernelVector(args{:})';

                    % Compute kernel fcn projection to H(m-1)
                    i.init(kexp.getKernelMatrix);
                    A = i.interpolate(Kbig');
                    F = atd.fxi(:,this.Used);
                    
%                     kexp.Ma = i.interpolate(atd.fxi(:,this.Used))';
%                     FunVis2D(kexp, atd);

                    % Cap too small norms!
                    phinormsq = max(abs(1 - sum(A.*Kbig,1)), this.PhiNormMin^2);
                    
                    % Squared versions
                    fDotPhiSq = (atd.fxi(:,free) - F*A).^2;
                    
                    if exp_mode 
                        this.boundOutput(kexp, fDotPhiSq, kexp_prev);
                        kexp_prev = kexp.clone;
                    end
                    
                    fDotPhiSq = sum(fDotPhiSq,1);
                    
                    % Compute relative error
                    [v, maxidx] = max(fDotPhiSq);
                    this.MaxErrors(cidx,siz) = sqrt(v);
                    this.relerr(cidx,siz) = sqrt(max(fDotPhiSq ./ atdfxinormsq(free)));
                    
                    % Check stopping conditions
                    if this.relerr(cidx,siz) <= this.MaxRelErr 
                        break;
                    end
                    
                    % For the orthonormal variant, divide by phinormsq
                    if ~this.UseOGA
                        fDotPhiSq = fDotPhiSq ./ phinormsq;
                        [v, maxidx] = max(fDotPhiSq);
                        this.Gain(cidx,siz) = v;
                    else
                        this.Gain(cidx,siz) = v / phinormsq(maxidx);
                    end
                    
                    this.extendExpansion(kexp, atd, free(maxidx));
                    siz = siz+1;
                    
                    off = this.HerrDecay(siz-1);
                    off2 = this.VKOGABound(siz-1);
                    
                    % @todo use cumsum at the end!
                    this.HerrDecay(cidx, siz) = off + this.Gain(cidx,siz);
                    % Only build improved bound if VKOGA is used
                    if ~this.UseOGA
                        this.VKOGABound(cidx, siz) = off2 + 1/max(phinormsq);
                    end
                    
                    %% Debug - change later to better stopping criteria not using the kernel expansion
                    if exp_mode
                        projpart_old = this.exp_output(Kbig, A, kexp, phinormsq, atd, ...
                            free, fDotPhiSq, maxidx, v, cidx, siz, this.Used, projpart_old, i, args);
                    end
                    
                    [~,id] = lastwarn;
                    if strcmp(id,'MATLAB:nearlySingularMatrix')
                        warncnt = warncnt + 1;
                        lastwarn('');
                    end
                end
                
                if warncnt == 2
                    siz = siz-2;
                end
                
                [~, siz] = min(this.relerr(cidx,1:siz));
                
                this.expsizes(cidx) = siz;
                
                if this.relerr(cidx,siz) < minerr
                    bestused = this.Used(1:siz);
                    minerr = this.relerr(cidx,siz);
                    bestcidx = cidx;
                end
                pi.step;
            end
            
            % Restore best configuration
            ec.applyConfiguration(bestcidx, kexp);
            ec.StateConfig.vBestConfigIndex = bestcidx;
            kexp.setCentersFromATD(atd, bestused);
            % Compute Ma
            i.init(kexp.getKernelMatrix);
            kexp.Ma = i.interpolate(atd.fxi(:,bestused))';
            this.Used = bestused;
            
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
                fprintf('VKOMP best kernel config index:%d (value:%e) for VKOMP with OGA=%s, exp-size:%d\n',...
                    bestcidx,g,og,size(bestMa,2));
            end
            
            pi.stop;
        end
    end
    
    methods(Access=private)
        function projpart = exp_output(this, Kbig, A, kexp, phinormsq, atd, ...
                free, fDotPhiSqAll, maxidx, v, cidx, siz, used, projpart_old, i, args)
            
            i.init(kexp.getKernelMatrix);
            kexp.Ma = i.interpolate(atd.fxi(:,this.Used))';
            afxi = kexp.evaluate(args{:});
            this.err(cidx,siz) = max(this.ErrorFun(afxi-atd.fxi(:,free)));
            projpart = 0;
            if KerMor.App.Verbose > 2
                projpart = -sum(sum(kexp.Ma .* atd.fxi(:,used),2));
                
                gain_direct = abs(projpart - projpart_old);
                gain_impl = this.Gain(cidx,siz);
                
                fprintf('Rel err:%e, gains dir: %e, impl: %e, rel. diff (to impl): %e, ||f^m||_H^2: %e\n',...
                    this.relerr(cidx,siz),gain_direct,gain_impl,(gain_direct-gain_impl)/gain_impl, kexp.NativeNorm^2);
                
                if KerMor.App.Verbose > 3
                    if isa(this.f,'function_handle')
                        FunVis2D(kexp, atd, [], this.f);
                    end
                    
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
                    atdfxi = atd.fxi.toMemoryMatrix;
                    plot(pxifree,fDotPhiSqAll,'r',pxi,atdfxi(yplotdim,:),'b',...
                        pxi,afxi(yplotdim,:),'g',pxifree(maxidx),v,'b.',cxi,cafxi,'g.');
                    title(sprintf('Extension choice %d',free(maxidx)));
                    lbl = '<f,phi>^2';
                    if this.UseOGA
                        lbl = '<f,phi~>^2';
                    end
                    legend(lbl,'f(x)','approx');
                    axis tight;
                    
                    subplot(2,3,5);
                    semilogy(this.Gain(cidx,1:siz));
                    title(sprintf('Gain from %d to %d basis functions: %e',...
                        siz-1,siz,this.Gain(cidx,siz)));
                    axis tight;
                    
                    subplot(2,3,6);
                    semilogy(this.HerrDecay(cidx,1:siz));
                    title(sprintf('Accumulated H projection error at %d basis functions',siz));
                    axis tight;
                    
                    pause;
                end
            end
        end
        
        function boundOutput(this, kexp, fDotPhiSqAll_comp, kexp_prev)
            if ~this.UseOGA && ~isempty(this.f) && KerMor.App.Verbose > 1 ...
                    && isa(this.f,'kernels.KernelExpansion')
                M = this.f.MBnd;
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
                hlp2 = this.f.NativeNorm^2 - this.HerrDecay(cidx, siz);
                pt.addRow(hlp1, hlp2);
                pt.display;
                if abs(hlp1-hlp2)/abs(hlp1) > 5e-4
                    warning('freak:out','Inconsistent ||f-f^m|| vs. ||f|| - accumulated gain');
                end
            end
        end
    end
end
