classdef MinMaxAdaptiveCWKA < approx.algorithms.BaseAdaptiveCWKA
% Adaptive component-wise kernel approximation algorithm
%
% @author Daniel Wirtz @date 2011-11-03
%
% See also: BaseApprox KernelApprox BaseAdaptiveCWKA
%
% @change{0,6,dw,2011-11-16} Fixed a bug that caused the best `M_\alpha` value to be overridden
% and thus messing up the kernel expansion. Now the best expansion for both gamma
% configurations and extensions is returned.
%
% @new{0,5,dw,2011-11-03} Added this class.

    properties(SetObservable)
        % Determines how many percent of the samples are checked as
        % possible center extensions. The samples are sorted descending by
        % the current error on each sample.
        %
        % Allowed range is `[0, 1]`. Zero means take only the maximum
        % error.
        %
        % @propclass{important} The size of the training data to use for
        % possible extension balances the approximation sharpness with the
        % computational costs.
        %
        % @default .1 @type double
        CheckMaxErrorPercent = .1;
    end
    
    properties(SetAccess=private)
        % Some debug data regarding the CheckMaxErrorPercent property
        %
        % First entry: Maximum % value of CheckMaxErrorPercent that would
        % have been effectively needed to obtain results (if =
        % CheckMaxErrorPercent, then it also could have been needed higher)
        % Second entry: Average over all percentage values
        CheckMaxErrorPercentData;
    end
    
    methods    
        function this = MinMaxAdaptiveCWKA
            this = this@approx.algorithms.BaseAdaptiveCWKA;
        end
                        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.MinMaxAdaptiveCWKA;
            
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
            
            numSamples = size(atd.xi,2);
            this.CheckMaxErrorPercentData = [];
            
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
                kexp.Ma = atd.fxi(:,inidx);
                used = inidx;
                
                %% Set kernel config - gammas
                this.setDistKernelConfig(kexp, dists(:,distidx));
                
                if KerMor.App.Verbose > 2
                    info = sprintf('Current \\beta_s: %f',kexp.Kernel.Gamma);
                    if this.pte
                        if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                            info = sprintf('%s, \\beta_t: %f',info,kexp.TimeKernel.Gamma);
                        end
                        if ~isa(kexp.ParamKernel,'kernels.NoKernel')
                            info = sprintf('%s, \\beta_\\mu: %f',info,kexp.ParamKernel.Gamma);
                        end
                    end
                end
                
                % Prepare new kernel matrix of size one
                K = data.MemoryKernelMatrix(1);
                K.BuildInverse = true;
                
                cnt = 2;
                % Start with minimum errror at one-center-expansion
                [sminerr, ~, minerrs] = this.getError(kexp, atd);
                this.MaxErrors(distidx,1) = sminerr;
                % Kernel expansion augmentation loop
                while true
                    % Compile the number of centers to try for extension
                    minerrs(used) = []; % Remove errors on used centers (centers are unique)
                    hlp = setdiff(1:numSamples,used);
                    [verb_out_only, sortidx] = sort(minerrs,'descend');
                    sortidx = hlp(sortidx);
                    l = length(sortidx);
                    if this.CheckMaxErrorPercent == 0
                        sortidx = sortidx(1);
                    else
                        sortidx = sortidx(1:ceil(l*this.CheckMaxErrorPercent));
                    end
                    
                    if KerMor.App.Verbose > 1
                        fprintf('Checking %d expansion extensions to size %d with minimum max error of %e\n',length(sortidx),cnt,verb_out_only(1))
                    end
                    
                    % Loop over all possible centers to add
                    localminerr = Inf;
                    for cidx = 1:length(sortidx)
                        sidx = sortidx(cidx);
                        % Set current center & get kernel vector
                        kexp.Centers.xi(:,cnt) = atd.xi(:,sidx);
                        if this.pte
                            if atd.hasTime
                                kexp.Centers.ti(cnt) = atd.ti(:,sidx);
                            end
                            if atd.hasParams
                                kexp.Centers.mui(:,cnt) = atd.mui(:,sidx);
                            end
                            PhiN = kexp.getKernelVector(atd.xi(:,sidx),atd.ti(sidx),atd.mui(:,sidx));
                        else
                            PhiN = kexp.getKernelVector(atd.xi(:,sidx));
                        end
                        Ktmp = K.clone;
                        Ktmp.augment(PhiN);

                        % Compute coefficients
                        this.CoeffComp.init(Ktmp, kexp);
                        this.computeCoeffs(kexp, atd.fxi(:,[used sidx]), [kexp.Ma zeros(size(atd.xi,1),1)]);
                        
                        % Get error
                        [val, ~, errs] = this.getError(kexp, atd);
                        if val < localminerr
                            minerr_sidx = sidx;
                            minerr_cidx = cidx;
                            
                            localminerr = val;
                            minerrs = errs;
                            
                            bestK = Ktmp;
                            innerbestMa = kexp.Ma;
                        end
                    end
                    if sminerr < localminerr && KerMor.App.Verbose > 1
                        fprintf('Warning. No considered center ext to size %d caused a decrease in error (prev: %e, current: %e)\n',cnt,sminerr, localminerr);
                    end
                    sminerr = localminerr;

                    K = bestK; % Use bestK as kernel matrix
                    % Assign best found center (given by minerr_sidx)
                    kexp.Centers.xi(:,cnt) = atd.xi(:,minerr_sidx);
                    if this.pte
                        if atd.hasTime
                            kexp.Centers.ti(cnt) = atd.ti(:,minerr_sidx);
                        end
                        if atd.hasParams
                            kexp.Centers.mui(:,cnt) = atd.mui(:,minerr_sidx);
                        end
                    end
                    % Set best coefficients from inner loop (needed to be set here as calls to
                    % compute coeffs could involve using current coefficients (SVR etc)
                    kexp.Ma = innerbestMa; 
                    
                    rel = sminerr / (norm(atd.fxi(minerr_sidx))+eps);
                    this.MaxErrors(distidx,cnt) = sminerr;
                    
                    this.CheckMaxErrorPercentData(end+1) = 100*minerr_cidx/numSamples;
                    if KerMor.App.Verbose > 1
                        fprintf('Minimum extended error %e achieved for center %d/%d (%2.2f%% of tr-data, real index:%d)\n',...
                            sminerr,minerr_cidx,length(sortidx),this.CheckMaxErrorPercentData(end),minerr_sidx);
                        if KerMor.App.Verbose > 2
                            doPlots;
                        end
                    end
                    
                    used(end+1) = minerr_sidx;%#ok

                    %% Stopping condition
                    if this.checkStop(cnt, rel, sminerr)
                        break;
                    end
                    cnt = cnt+1;
                end
                
                % Store approx if this gamma config gave a better solution
                if sminerr < minerr
                    bestMa = kexp.Ma;
                    bestdistidx = distidx;
                    bestused = used;
                    minerr = sminerr;
                end
            end
            
            %% Set kernel expansion config of best found one
            kexp.Ma = bestMa;
            this.setDistKernelConfig(kexp, dists(:,bestdistidx));
            kexp.Centers.xi = atd.xi(:,bestused);
            if this.pte
                if atd.hasTime
                    kexp.Centers.ti = atd.ti(bestused);
                end
                if atd.hasParams
                    kexp.Centers.mui = atd.mui(:,bestused);
                end
            end
            
            % Statistics
            this.CheckMaxErrorPercentData = [max(this.CheckMaxErrorPercentData) mean(this.CheckMaxErrorPercentData)];
                
            if KerMor.App.Verbose > 1
                figure;
                semilogy(1:size(this.MaxErrors,2),this.MaxErrors);
            end
        
            function doPlots
                figure(1);
                fprintf('Exp-size: %d, max-err: %5.15f, rel-err: %5.15f, %s\n', cnt, sminerr, rel,info);
                pos = [1 3];
                if KerMor.App.Verbose > 3
                    pos = 1;
                end
                subplot(2,2,pos);
                plot(1:length(errs),errs,'r',used,val,'r*',minerr_sidx,val,'b*');
                title(info);
                subplot(2,2,pos+1);
                hold off;
                plot([atd.Box.xmin, atd.Box.xmax],'black');
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
                    plot([atd.Box.tmin, atd.Box.tmax],'black');
                    axis tight;
                    hold on;
                    if size(kexp.Centers.ti,2) > 1
                        plot(kexp.Centers.ti(:,1:end-1),'.','MarkerSize',2);
                    end
                    plot(kexp.Centers.ti(:,end),'r*','MarkerSize',3);
                    subplot(2,2,4); hold off;
                    plot([atd.Box.mumin, atd.Box.mumax],'black');
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


