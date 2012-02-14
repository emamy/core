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
        Herr;
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
            [~, initialidx] = max(sum(atd.fxi,1));
            % Supremum
            %[v, idx] = max(max(abs(atd.fxi),[],1));
            
            i = general.interpolation.KernelInterpol;
            total = 1:size(atd.xi,2);
            
            d = this.getDists(atd);
            this.err = zeros(this.NumGammas,this.MaxExpansionSize);
            this.Herr = this.err;
            this.Herr(:,1) = sqrt(sum(atd.fxi(:,initialidx).^2));
            
            minHerr = Inf;
            for gidx = 1:size(d,2)
                
                this.setDistKernelConfig(kexp, d(:,gidx))
                
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
                i.init(kexp.getKernelMatrix);
                kexp.Ma = i.interpolate(atd.fxi(:,used))';
                afxi = kexp.evaluate(atd.xi);
                this.err(gidx,1) = max(sqrt(sum((afxi-atd.fxi).^2)));
                
                siz = 1;
                while siz < this.MaxExpansionSize
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
                    phinormsq = 1 - 2*sum(A.*Kbig,1) + sum(A.*(K*A),1);
                    %A = A./repmat(phinorm,siz,1);
                    %semilogy(phinorm);

                    F = atd.fxi(:,used);
                    fDotPhiSqAll = (atd.fxi(:,free) - F*A).^2;

                    % Sum
                    fDotPhiSqAll = sum(fDotPhiSqAll,1);
                    % Supremum
                    %fDotPhiSqAll = max(abs(fDotPhiSqAll),[],1);

    %                 [~, maxidxphi] = max(sum(fDotPhiSqAll,1)./phinorm);
    %                 fprintf('Max idx with/without norming: %d/%d\n',maxidx,maxidxphi);
    %                 maxidx = maxidxphi;

                    [v, maxidx] = max(fDotPhiSqAll .* phinormsq);
                    %[v, maxidx] = max(fDotPhiSqAll);
                    
                    %plot(fDotPhiSqAll)
                    this.extendExpansion(kexp, atd, free(maxidx));
                    used(end+1) = free(maxidx);%#ok

                    siz = siz+1;
                    
                    this.Herr(gidx,siz) = this.Herr(gidx,siz-1)-v*phinormsq(maxidx);
                    
                    %% Debug
                    i.init(kexp.getKernelMatrix);
                    kexp.Ma = i.interpolate(atd.fxi(:,used))';
                    afxi = kexp.evaluate(atd.xi);
                    this.err(gidx,siz) = max(sqrt(sum((afxi-atd.fxi).^2)));
                end
                
                if this.err(gidx,end) < minHerr
                    bestused = used;
                    bestMa = kexp.Ma;
                    bestDist = d(:,gidx);
                    minHerr = this.err(gidx,end);
                    bestgidx = gidx;
                end
                %i.init(kexp.getKernelMatrix);
                %kexp.Ma = i.interpolate(atd.fxi(:,used))';
            end
            bestgidx
            this.setDistKernelConfig(kexp, bestDist);
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