classdef AdaptiveCompWiseKernelApprox < approx.algorithms.BaseKernelApproxAlgorithm
% Adaptive component-wise kernel approximation algorithm
%
% @author Daniel Wirtz @date 2011-03-31
%
% See also: BaseApprox KernelApprox
%
% @change{0,5,dw,2011-09-09} Fixed setters for MaxRelErr and
% MaxAbsErrFactor, now setting values to 'this' instead of 'kexp'
%
% @change{0,5,dw,2011-07-28} Changed the algorithm part so that it can also work with
% kernels.KernelExpansion instead of only on kernels.ParamTimeKernelExpansion.
%
% @new{0,5,dw,2011-07-07} Moved the old approx.AdaptiveCompWiseKernelApprox class to this class.
%
% @change{0,4,dw,2011-05-31} Added new experimental properties @ref MinGFactor and @ref MaxGFactor.
%
% @change{0,4,dw,2011-05-19} Disconnected the Approx classes from taking a BaseModel instance at
% approx computation. This way external tools can use the approximation algorithms, too.
%
% @change{0,3,sa,2011-04-21} Implemented Setters for all the properties
% other than NumGammas and ValidationPercent
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @change{0,3,dw,2011-04-14}
% - Implemented some setters
% - New property approx.algorithms.AdaptiveCompWiseKernelApprox.ValidationPercent enabling a validation set
% to check for best gammas
%
% @change{0,3,dw,2011-04-06} Now works with models that dont have any
% parameters.
%
% @new{0,3,dw,2011-04-01} Added this class.
    
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
        
        % Percentage `p` of the training data to use as validation data
        %
        % Admissible values are `p\in]0,\frac{1}{2}]`.
        %
        % @propclass{optional} 
        %
        % @default .2
        ValidationPercent = .2;
        
        % Value for initial Gamma choice.
        %
        % @propclass{experimental} 
        gameps = .6;
        
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
        MinGFactor = .05;
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
        function this = AdaptiveCompWiseKernelApprox
            this = this@approx.algorithms.BaseKernelApproxAlgorithm;
            
            % Register default property changed listeners
            this.registerProps('MaxExpansionSize','NumGammas','ValidationPercent',...
                'gameps','MaxRelErr','MaxAbsErrFactor','ErrFun','MaxGFactor','MinGFactor');
        end
                        
%         function target = clone(this)
%             % Clones the instance.
%             
%             % Create instance as this is the final class so far. If
%             % subclassed, this clone method has to be given an additional
%             % target argument.
%             target = approx.algorithms.AdaptiveCompWiseKernelApprox;
%             
%             target = clone@approx.KernelApprox(this, target);
%             
%             %this.cloneLocalProps(target,mfilename('class'));
%             % copy local props
%             copy.MaxExpansionSize = kexp.MaxExpansionSize;
%             copy.NumGammas = this.NumGammas;
%             copy.gameps = this.gameps;
%             copy.MaxRelErr = kexp.MaxRelErr;
%             copy.MaxAbsErrFactor = kexp.MaxAbsErrFactor;
%             copy.MaxErrors = kexp.MaxErrors;
%             copy.ValidationPercent = this.ValidationPercent;
%             copy.ErrFun = this.ErrFun;
%             copy.effabs = this.effabs;
%         end
    end
    
    methods(Access=protected, Sealed)
        function detailedComputeApproximation(this, kexp, xi, ti, mui, fxi)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
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
            
            % Extract validation data
            vnum = round(size(xi,2)*this.ValidationPercent);
            sel = round(linspace(1,size(xi,2),vnum));
            vdxi = xi(:,sel);
            xi(:,sel) = [];
            vdfxi = fxi(:,sel);
            fxi(:,sel) = [];
            
            if pte
                vdti = ti(sel);
                ti(sel) = [];
            end
            
            % check if training data contains parameters
            nx = general.NNTracker;
            if pte
                nt = general.NNTracker;
            end
            
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
                    vdmui = mui(:,sel);
                    mui(:,sel) = [];

                    np = general.NNTracker;

                    [BPmin, BPmax] = general.Utils.getBoundingBox(mui);
                    bpdia = norm(BPmax - BPmin);

                    thecenter = [thecenter; (BPmin + BPmax)/2];
                    B = [B; mui];
                else
                    hasparams = false;
                    vdmui = [];
                end
            end
            
            %% Select initial center x0
            % Strategy: Take the point that is closest to the bounding
            % box center!
            A = repmat(thecenter, 1, size(xi,2));
            [dummy, inIdx] = min(sum((A-B).^2,1));
            clear A B;
            
            kexp.Centers.xi = xi(:,inIdx);
            nx.addPoint(xi(:,inIdx)); % Add points to nearest neighbor trackers (for gamma comp)
            if pte
                kexp.Centers.ti = ti(inIdx);
                nt.addPoint(ti(inIdx));
                if hasparams
                    kexp.Centers.mui = mui(:,inIdx);
                    np.addPoint(kexp.Centers.mui);
                else
                    kexp.Centers.mui = [];
                end
            end
            
            %% Set up initial expansion
            used = inIdx;
            kexp.Ma = fxi(:,inIdx);
            
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
            
            % Init for the loop
            minerr = Inf;             
            gt = []; gp = [];
            for idx = 1:size(dists,2)
                gx = kexp.Kernel.setGammaForDistance(dists(1,idx),this.gameps);
                if pte
                    if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                        gt = kexp.TimeKernel.setGammaForDistance(dists(2,idx),this.gameps);
                    end
                    if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                        gp = kexp.ParamKernel.setGammaForDistance(dists(3,idx),this.gameps);
                    end
                end
                
                if pte
                    val = errfun(fxi,kexp.evaluate(xi,ti,mui));
                else
                    val = errfun(fxi,kexp.evaluate(xi));
                end
                if val < minerr
                    minerr = val;
                    bestgx = gx;
                    bestgt = gt;
                    bestgp = gp;
                    bestdistidx = idx;
                end
            end
            %% Assign best values
            kexp.Kernel.Gamma = bestgx;
            if pte
                if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                    kexp.TimeKernel.Gamma = bestgt;
                end
                if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                    kexp.ParamKernel.Gamma = bestgp;
                end
            end
            if KerMor.App.Verbose > 1
                fprintf('Initial gammas: SK:%e, TK:%e, PK:%e\n',bestgx,bestgt,bestgp);
            end
            
            %% Outer control loop
            cnt = 1;
            % Stopping condition preps
            this.effabs = this.MaxAbsErrFactor * max(abs(fxi(:)));
            
            % Keep track of maximum errors
            this.MaxErrors = zeros(1,this.MaxExpansionSize);
            exception = false;
            while ~exception
                
                %% Determine maximum error over training data
                if pte
                    fhat = kexp.evaluate(xi,ti,mui);
                else
                    fhat = kexp.evaluate(xi);
                end
                [val, maxidx, errs] = errfun(fxi,fhat);
                rel = val / (norm(fxi(maxidx))+eps);
                this.MaxErrors(cnt) = val;
                
                %% Verbose stuff
                if KerMor.App.Verbose > 2
                    doPlots;
                end
                
                %% Stopping condition
                if this.checkStop(cnt, rel, val)
                    break;
                end
                
                % Add maxidx to list of used centers
                used(end+1) = maxidx;%#ok
                
                %% Extend centers
                kexp.Centers.xi(:,end+1) = xi(:,maxidx);
                if pte
                    kexp.Centers.ti(end+1) = ti(maxidx);
                    if hasparams
                        kexp.Centers.mui(:,end+1) = mui(:,maxidx);
                        np.addPoint(mui(:,maxidx));
                    end
                    nt.addPoint(ti(maxidx));
                end
                % Add points to nearest neighbor trackers (for gamma comp)
                nx.addPoint(xi(:,maxidx));
                
                %% Compute new approximation
                olddists = dists;
                xdists = sort([dfun(nx.getMinNN, this.MaxGFactor*bxdia) olddists(1,bestdistidx)]);
                dists = xdists;
                if pte
                    tdists = sort([dfun(nt.getMinNN, this.MaxGFactor*btdia) olddists(2,bestdistidx)]);
                    dists = [dists; tdists];%#ok
                    if hasparams
                        minnn = bpdia/this.NumGammas;
                        if ~isinf(np.getMinNN)
                            minnn = np.getMinNN;
                        end
                        pdists = sort([dfun(minnn, this.MaxGFactor*bpdia) olddists(3,bestdistidx)]);
                        dists = [dists; pdists];%#ok
                    end
                end
                
                if KerMor.App.Verbose > 2
                    %dists%#ok
                end
                
                minerr = Inf; gt = []; gp = [];
                for gidx = 1:size(dists,2)
                    d = dists(:,gidx);
                    % Update Kernels Gamma values
                    gx = kexp.Kernel.setGammaForDistance(d(1),this.gameps);
                    if KerMor.App.Verbose > 2
                        %fprintf('Kernels - Sys:%10f => gamma=%f',d(1),gx);
                        fprintf('xg:%.5e',gx);
                    end
                    if pte
                        if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                            gt = kexp.TimeKernel.setGammaForDistance(d(2),this.gameps);
                            if KerMor.App.Verbose > 2
                                %fprintf(', Time:%10f => gamma=%10f',d(2),gt);
                                fprintf(', tg:%.5e',gx);
                            end
                        end
                        if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                            gp = kexp.ParamKernel.setGammaForDistance(d(3),this.gameps);
                            if KerMor.App.Verbose > 2
                                fprintf(', pg=%.5e',gp);
                            end
                        end
                    end
                    
                    %% Compute coefficients
                    %warning('off','MATLAB:nearlySingularMatrix');
                    % Call coeffcomp preparation method and pass kernel matrix
                    K = data.MemoryKernelMatrix(kexp.getKernelMatrix);
                    %K.UseLU = true;
                    %K.BuildInverse = true;
                    this.CoeffComp.init(K);
                    
                    % Call protected method
                    try
                        this.computeCoeffs(kexp, fxi(:,used), [kexp.Ma zeros(size(fxi,1),1)]);
                    catch ME
                        if strcmp(ME.identifier,'KerMor:coeffcomp:failed')                            
                            exception = true;
                            ex = ME;
                            break;
                        else
                            rethrow(ME);
                        end
                    end
                    %warning('on','MATLAB:nearlySingularMatrix');
                    
                    % get error on training data
                    if pte
                        fhat = kexp.evaluate(xi,ti,mui);
                    else
                        fhat = kexp.evaluate(xi);
                    end
                    val = errfun(fxi, fhat);
                    impro = (val / minerr) * 100;
                    
                    % get error on validation set
                    if pte
                        fvali = kexp.evaluate(vdxi,vdti,vdmui);
                    else
                        fvali = kexp.evaluate(vdxi);
                    end
                    vdval = errfun(vdfxi, fvali);
                    if val < minerr
                        minerr = val;
                        bestgx = gx;
                        bestgt = gt;
                        bestgp = gp;
                        bestMa = kexp.Ma;
                        bestdistidx = idx;
                        if KerMor.App.Verbose > 2
                            fprintf(' b: %.5e, %3.2f%%',val,impro);
                        end
                    else
                        if KerMor.App.Verbose > 2
                            fprintf(' w: %.5e, %3.2f%%',val,impro);
                        end
                    end
                    
                    if KerMor.App.Verbose > 2
                        fprintf(' ||Ma||:%.5e, vd-err:%.5e, prod:%.5e\n',sum(sqrt(sum(kexp.Ma.^2,1))),vdval,val*vdval);
                    end
                end
                if KerMor.App.Verbose > 2
                    fprintf('\n');
                end
                
                %% Assign best values
                kexp.Kernel.Gamma = bestgx;
                if pte
                    if ~isa(kexp.TimeKernel,'kernels.NoKernel')
                        kexp.TimeKernel.Gamma = bestgt;
                    end
                    if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
                        kexp.ParamKernel.Gamma = bestgp;
                    end
                end
                kexp.Ma = bestMa;
                
                if KerMor.App.Verbose > 1
                    fprintf('-- It: %d ---- Minerr: %f ----- Best values: System:%f, Time:%f, Param:%f ----------\n',cnt,minerr,bestgx,bestgt,bestgp);
                end
                
                cnt = cnt+1;
            end
            
            if exception
                fprintf('Adaptive approximation generation stopped due to exception.\n')
                cprintf('red',ex.getReport);
            end
            
            if KerMor.App.Verbose > 1
                figure;
                plot(this.MaxErrors,'r');
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
                fprintf('Max error over training data: %5.20f (relative: %10.20f)\n',val,rel);
                pos = [1 3];
                if KerMor.App.Verbose > 3
                    pos = 1;
                end
                subplot(2,2,pos);
                plot(1:length(errs),errs,'r',used,val,'r*',maxidx,val,'b*');
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
        function set.ValidationPercent(this, value)
            if ~isposrealscalar(value) || value > .5
                error('The value must be a positive scalar inside the interval ]0,.5[');
            end
            this.ValidationPercent = value;
        end
        
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
                              
        function set.gameps(this, value)
            if ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.gameps = value;
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
    
    methods(Access=private)
        function bool = checkStop(this, cnt, rel, val)
            % Checks the stopping conditions 
            bool = false;
            if cnt == this.MaxExpansionSize
                fprintf('AdaptiveCompWiseKernelApprox finished. Max expansion size %d reached.\n',this.MaxExpansionSize);
                bool = true;
            elseif rel < this.MaxRelErr
                fprintf('AdaptiveCompWiseKernelApprox finished. Relative error %.7e < %.7e\n',rel,this.MaxRelErr);
                bool = true;
            elseif val < this.effabs
                fprintf('AdaptiveCompWiseKernelApprox finished. Absolute error %.7e < %.7e\n',val,this.effabs);
                bool = true;
            end
        end
    end
end


