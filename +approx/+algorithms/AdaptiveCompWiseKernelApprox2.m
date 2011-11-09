classdef AdaptiveCompWiseKernelApprox2 < approx.algorithms.BaseAdaptiveCWKA
% Adaptive component-wise kernel approximation algorithm, improved version
%
% This adaptive algorithm loops over different kernel configurations as the
% main outer loop instead of changing kernel spaces inside the extension
% loop. The tried gamma values are fixed and are not adopted to the current
% kernel expansion (range-wise).
%
% @author Daniel Wirtz @date 2011-09-09
%
% See also: BaseApprox KernelApprox
%
% @change{0,5,dw,2011-11-02}
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
% - Moved common properties for adaptive algorithms to new base class BaseAdaptiveCWKA and
% using it's provided methods in order to have a more compact algorithm representation.
%
% @new{0,5,dw,2011-09-09} Added this class.
%
% @change{0,5,dw,2011-11-03} Using BaseAdaptiveCWKA interface now.

    methods    
        function this = AdaptiveCompWiseKernelApprox2
            this = this@approx.algorithms.BaseAdaptiveCWKA;
        end

        function copy = clone(this)
            % Clones the instance.
            copy = approx.algorithms.AdaptiveCompWiseKernelApprox2;
            copy = clone@approx.algorithms.BaseAdaptiveCWKA(this, copy);
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedAdaptiveApproximation(this, kexp, atd)
            % Performs adaptive approximation generation.
            %
            % Parameters:
            % kexp: The kernel expansion. @type kernels.KernelExpansion
            % atd: The approximation training data instance @type data.ApproxTrainData
            
            [c, inidx] = this.getInitialCenter(atd);
            
            % Choose gammas
            dists = this.getDists(atd);
            distnum = size(dists,2);
            
            %% Init for the loop
            % Keep track of maximum errors
            this.MaxErrors = zeros(distnum,this.MaxExpansionSize);
            
            minerr = Inf;        
            for distidx = 1:distnum
                %% Set up initial expansion
                kexp.Centers.xi = c(1:size(atd.xi,1));
                if this.pte
                    if atd.hasTime
                        kexp.Centers.ti = c(atd.tOff);
                    end
                    if atd.hasParams
                        kexp.Centers.mui = c(atd.muOff:end);
                    end
                end
                kexp.Ma = atd.fxi(:,inIdx);
                bestMa = kexp.Ma;
                used = inidx;
                
                %% Set kernel config - gammas
                this.setDistKernelConfig(kexp, dists(:,distidx));
                
                if KerMor.App.Verbose > 2
                    info = sprintf('Current \\beta_s: %f',kexp.Kernel.Gamma);
                    if pte
                        if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                            info = sprintf('%s, \\beta_t: %f',info,kexp.TimeKernel.Gamma);
                        end
                        if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                            info = sprintf('%s, \\beta_\\mu: %f',info,kexp.ParamKernel.Gamma);
                        end
                    end
                end
                
                %% Adaptive extension loop
                % prepare new kernel matrix of size one
                K = data.MemoryKernelMatrix(1);
                K.BuildInverse = true;
                %K.UseLU = true;
                this.CoeffComp.init(K);
                
                cnt = 1;
                while true
                    %% Determine maximum error over training data
                    [val, maxidx, errs] = this.getError(kexp, atd);
                    % |@todo crap computation
                    rel = val / (norm(atd.fxi(maxidx))+eps);
                    this.MaxErrors(distidx,cnt) = val;
                    
                    if KerMor.App.Verbose > 2
                        doPlots;
                    end

                    %% Stopping condition
                    if this.checkStop(cnt, rel, val)
                        break;
                    end
                    
                    %% Extend centers
                    used(end+1) = maxidx;%#ok
                    this.extendExpansion(kexp, atd, maxidx);
                    if this.pte
                        PhiN = kexp.getKernelVector(atd.xi(:,maxidx),atd.ti(maxidx),atd.mui(:,maxidx));
                    else
                        PhiN = kexp.getKernelVector(atd.xi(:,maxidx));
                    end
                    K.augment(PhiN);
                    cnt = cnt+1;
                                        
                    %% Compute coefficients
                    this.computeCoeffs(kexp, atd.fxi(:,used), [kexp.Ma zeros(d,1)]);
                end
                
                %% Store indices if better approx is found
                if val < minerr
                    minerr = val;
                    bestMa = kexp.Ma;
                    bestdistidx = distidx;
                end
                
            end
            
            %% Set kernel expansion config of best found one
            kexp.Ma = bestMa;
            this.setDistKernelConfig(kexp, dists(:,bestdistidx));
                
            if KerMor.App.Verbose > 1
                figure;
                semilogy(1:size(this.MaxErrors,2),this.MaxErrors);
            end
        
            function doPlots
                figure(1);
                fprintf('Exp-size: %d, max-err: %5.15f, rel-err: %5.15f, %s\n',cnt, val, rel,info);
                pos = [1 3];
                if KerMor.App.Verbose > 3
                    pos = 1;
                end
                subplot(2,2,pos);
                plot(1:length(errs),errs,'r',used,val,'r*',maxidx,val,'b*');
                title(info);
                subplot(2,2,pos+1);
                hold off;
                plot([BXmin, BXmax],'black');
                axis tight;
                hold on;
                %plot((BXmin+BXmax)/2,'g');
                if size(kexp.Centers.xi,2) > 1
                    plot(kexp.Centers.xi(:,1:end-1),'.','MarkerSize',2);
                end
                plot(kexp.Centers.xi(:,end),'r*','MarkerSize',3);
                
                % Plot params & time also
                if KerMor.App.Verbose > 3
                    subplot(2,2,3); hold off;
                    plot([Btmin, Btmax],'black');
                    axis tight;
                    hold on;
                    if size(kexp.Centers.ti,2) > 1
                        plot(kexp.Centers.ti(:,1:end-1),'.','MarkerSize',2);
                    end
                    plot(kexp.Centers.ti(:,end),'r*','MarkerSize',3);
                    subplot(2,2,4); hold off;
                    plot([BPmin, BPmax],'black');
                    axis tight;
                    hold on;
                    if size(kexp.Centers.mui,2) > 1
                        plot(kexp.Centers.mui(:,1:end-1),'.','MarkerSize',2);
                    end                  
                    pause;
                end
            end
        end 
    end
end


