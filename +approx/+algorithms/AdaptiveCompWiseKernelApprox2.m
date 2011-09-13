classdef AdaptiveCompWiseKernelApprox2 < approx.algorithms.BaseKernelApproxAlgorithm
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
% @new{0,5,dw,2011-09-09} Added this class.
    
    properties(SetObservable)
        % The maximum size of the expansion to produce.
        %
        % Equals the maximum number of iterations to perform during
        % adaptive approximation computation as each iteration yields a new
        % center.
        %
        % @propclass{alglimit} 
        % Some text describing the importance of this property.
        %
        % @default 200
        MaxExpansionSize = 200;
        
        % The number of different Gamma values to try.
        %
        % @propclass{important} 
        NumGammas = 10;
        
        % Stopping condition property. Maximum relative error that may occur
        %
        % @propclass{critical}
        %
        % @default 1e-5
        MaxRelErr = 1e-5;
        
        % Stopping condition property. 
        %
        % Factor for maximum absolute error that may occur. Factor means that this value will be
        % multiplied by the maximum value of the approximation training data f-Values used to train
        % this approximation. This way, one defines the fraction of the maximum error that is
        % maximally allowed.
        %
        % @propclass{critical}
        %
        % @default 1e-5
        MaxAbsErrFactor = 1e-5;
        
        % The error functional to use
        % 1 = `L^\infty`-Error (max diff over all vecs & dimensions)
        % 2 = `L^2`-Error (vector-wise L^2, then max)
        %
        % @propclass{experimental} 
        ErrFun = 1;
        
        % 'dfun(this.MinGFactor*bxdia, kexp.MaxGFactor*bxdia);'
        % @propclass{experimental}
        MaxGFactor = 1;
        
        % 'dfun(this.MinGFactor*bxdia, kexp.MaxGFactor*bxdia);'
        % @propclass{experimental}
        MinGFactor = .1;
    end
    
    properties(Transient, SetAccess=private)
        % Contains the maximum errors for each iteration/center extension step performed by the last
        % run of this algorithm.
        MaxErrors = [];
    end
    
    properties(Transient, Access=private)
        effabs;
    end
    
    methods    
        function this = AdaptiveCompWiseKernelApprox2
            this = this@approx.algorithms.BaseKernelApproxAlgorithm;
            
            % Register default property changed listeners
            this.registerProps('MaxExpansionSize','NumGammas',...
                'MaxRelErr','MaxAbsErrFactor','ErrFun','MaxGFactor','MinGFactor');
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedComputeApproximation(this, kexp, xi, ti, mui, fxi)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
            gammaeps = .1;
            d = size(fxi,1);
            
            dfun = @logsp; % gamma distances comp fun (linsp / logsp)
            if this.ErrFun == 1
                errfun = @getLInftyErr; % L^inf error function
            else
                errfun = @getL2Err; % L^2 error function
            end
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            if ~isa(kexp, 'kernels.KernelExpansion')
                error('Approximation method works only for kernel expansions.');
            elseif ~isa(kexp.Kernel,'kernels.GaussKernel')
                error('The state kernel has to be a Gaussian for this approximation algorithm so far');
            end
            pte = isa(kexp,'kernels.ParamTimeKernelExpansion');
            if pte && ((~isa(kexp.TimeKernel,'kernels.GaussKernel') && ~isa(kexp.TimeKernel,'kernels.NoKernel')) || ...
                    (~isa(kexp.ParamKernel,'kernels.GaussKernel') && ~isa(kexp.ParamKernel,'kernels.NoKernel')))
                error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
            end
            
            %% Initializations
            % Stopping condition preps
            this.effabs = this.MaxAbsErrFactor * max(abs(fxi(:)));
            
            %% Compute bounding boxes & center
            [BXmin, BXmax] = general.Utils.getBoundingBox(xi);
            thecenter = (BXmin+BXmax)/2;
            B = xi;
            % Get bounding box diameters
            bxdia = norm(BXmax - BXmin);
            if pte
                Btmin = min(ti); Btmax = max(ti);
                btdia = Btmax-Btmin;
                thecenter = [thecenter; (Btmin+Btmax)/2];
                B = [B; ti];
                
                %% Check if params are used
                if ~isempty(mui)
                    hasparams = true;

                    [BPmin, BPmax] = general.Utils.getBoundingBox(mui);
                    bpdia = norm(BPmax - BPmin);

                    thecenter = [thecenter; (BPmin + BPmax)/2];
                    B = [B; mui];
                else
                    hasparams = false;
                end
            end
            
            %% Select initial center x0
            % Strategy: Take the point that is closest to the bounding
            % box center!
            A = repmat(thecenter, 1, size(xi,2));
            [dummy, inIdx] = min(sum((A-B).^2,1));
            clear A B;
            
            %% Choose initial gammas
            xdists = dfun(this.MinGFactor*bxdia, this.MaxGFactor*bxdia);
            dists = xdists;
            if pte
                tdists = dfun(this.MinGFactor*btdia, this.MaxGFactor*btdia);
                dists = [dists; tdists];
                if hasparams
                    pdists = dfun(this.MinGFactor*bpdia, this.MaxGFactor*bpdia);
                    dists = [dists; pdists];
                end
            end
            distnum = size(dists,2);
            
            %% Init for the loop
            % Keep track of maximum errors
            this.MaxErrors = zeros(distnum,this.MaxExpansionSize);
            
            minerr = Inf;             
            for distidx = 1:distnum
                %% Set up initial expansion
                kexp.Centers.xi = xi(:,inIdx);
                if pte
                    kexp.Centers.ti = ti(inIdx);
                    if hasparams
                        kexp.Centers.mui = mui(:,inIdx);
                    else
                        kexp.Centers.mui = [];
                    end
                end
                kexp.Ma = fxi(:,inIdx);
                used = inIdx;
                
                %% Set kernel config - gammas
                kexp.Kernel.setGammaForDistance(dists(1,distidx),gammaeps);
                if pte
                    if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                        kexp.TimeKernel.setGammaForDistance(dists(2,distidx),gammaeps);
                    end
                    if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                        kexp.ParamKernel.setGammaForDistance(dists(3,distidx),gammaeps);
                    end
                end
                
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
                exception = false;
                while ~exception
                    
                    
                    %% Determine maximum error over training data
                    if pte
                        fhat = kexp.evaluate(xi,ti,mui);
                    else
                        fhat = kexp.evaluate(xi);
                    end
                    [val, maxidx, errs] = errfun(fxi,fhat);
                    % |@todo crap computation
                    rel = val / (norm(fxi(maxidx))+eps);
                    this.MaxErrors(distidx,cnt) = val;
                    
%                     %% Verbose stuff
                    if KerMor.App.Verbose > 2
                        doPlots;
                    end

                    %% Stopping condition
                    if val < this.effabs || rel < this.MaxRelErr || cnt == this.MaxExpansionSize
                        break;
                    end
                    
                    %% Extend centers
                    used(end+1) = maxidx;%#ok
                    kexp.Centers.xi(:,end+1) = xi(:,maxidx);
                    if pte
                        kexp.Centers.ti(end+1) = ti(maxidx);
                        if hasparams
                            kexp.Centers.mui(:,end+1) = mui(:,maxidx);
                        end
                    end
                    if pte
                        PhiN = kexp.getKernelVector(xi(:,maxidx),ti(maxidx),mui(:,maxidx));
                    else
                        PhiN = kexp.getKernelVector(xi(:,maxidx));
                    end
                    K.augment(PhiN);
                    cnt = cnt+1;
                    
%                     this.CoeffComp.Vis = 3 * cnt > 6;
%                     if cnt > 8
%                         this.CoeffComp.Vis = 3;
%                     end
                    
                    %% Compute coefficients
                    %warning('off','MATLAB:nearlySingularMatrix');
                    % Call protected method
                    try
                        this.computeCoeffs(kexp, fxi(:,used), [kexp.Ma zeros(d,1)]);
                    catch ME
                        if strcmp(ME.identifier,'KerMor:coeffcomp:failed')                            
                            exception = true;
                            ex = ME;
                            break;
                        else
                            rethrow(ME);
                        end
                    end
                end
                
                %% Store indices if better approx is found
                if val < minerr || isempty(bestMa)
                    minerr = val;
                    bestMa = kexp.Ma;
                    bestdistidx = distidx;
                end
                
            end
            
            if exception
                fprintf('Adaptive approximation generation stopped due to exception.\n')
                cprintf('red',ex.getReport);
            end
            
            %% Set kernel expansion config of best found one
            kexp.Ma = bestMa;
            kexp.Kernel.setGammaForDistance(dists(1,bestdistidx),gammaeps);
            if pte
                if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                    kexp.TimeKernel.setGammaForDistance(dists(2,bestdistidx),gammaeps);
                end
                if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                    kexp.ParamKernel.setGammaForDistance(dists(3,bestdistidx),gammaeps);
                end
            end
                
            if KerMor.App.Verbose > 1
                figure;
                semilogy(1:size(this.MaxErrors,2),this.MaxErrors);
            end
        
            function [val,idx,errs] = getLInftyErr(a,b)
                % computes the 'L^\infty'-approximation error over the
                % training set for the current approximation
                
                %errs = max(abs((a-b) ./ (a+eps)));
                errs = max(abs(a-b),[],1);
                [val, idx] = max(errs);
            end
            
            function [val,idx,errs] = getL2Err(a,b)
                % computes the 'L^\infty'-approximation error over the
                % training set for the current approximation
                errs = sqrt(sum((a-b).^2,1));
                [val, idx] = max(errs);
            end
            
            function linsp(from, to)%#ok
                d = linspace(from,to,this.NumGammas);
            end
            
            function d = logsp(from, to)
                d = logspace(log10(from),log10(to),this.NumGammas);
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
    
    %% Getter & Setter
    methods
        function set.NumGammas(this,value)
            if ~isposintscalar(value)
                error('Value must be a positive integer.');
            end
            this.NumGammas = value;
        end
        
        function set.MaxExpansionSize(this, value)
            if ~isposintscalar(value)
                error('Value must be a positive integer.');
            end
            this.MaxExpansionSize = value;
        end
        
        function set.MaxRelErr(this, value)
            if ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.MaxRelErr = value;
        end
        
        function set.MaxAbsErrFactor(this, value)
            if ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.MaxAbsErrFactor = value;
        end
        
        function set.ErrFun(this, value)
            if ~(value == 1 || value == 2)
                error('Value must be either integer 1 or 2.');
            end
            this.ErrFun = value;
        end
    end
    
end


