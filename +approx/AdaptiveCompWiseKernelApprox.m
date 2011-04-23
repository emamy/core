classdef AdaptiveCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
    % Adaptive component-wise kernel approximation algorithm
    %
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % See also: BaseApprox BaseCompWiseKernelApprox
    %
    % @change{0,3,sa,2011-04-21} Implemented Setters for all the properties
    % other than NumGammas and ValidationPercent
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @change{0,3,dw,2011-04-14}
    % - Implemented some setters
    % - New property approx.AdaptiveCompWiseKernelApprox.ValidationPercent enabling a validation set
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
        NumGammas = 20;
        
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
            this = this@approx.BaseCompWiseKernelApprox;
            
            % Register default property changed listeners
            this.registerProps('MaxExpansionSize','NumGammas','ValidationPercent',...
                'gameps','MaxRelErr','MaxAbsErrFactor','ErrFun');
        end
        
        function approximateCoreFun(this, model)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
            %% Experimental settings
            fac = 50;
            minfac = .01; % min factor for BB diameters at initial gamma choice
            dfun = @logsp; % gamma distances comp fun (linsp / logsp)
            if this.ErrFun == 1
                errfun = @getLInftyErr; % L^inf error function
            else
                errfun = @getL2Err; % L^2 error function
            end
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            if ~isa(this.SystemKernel,'kernels.GaussKernel') || ...
                    (~isa(this.TimeKernel,'kernels.GaussKernel') && ~isa(this.TimeKernel,'kernels.NoKernel')) || ...
                    (~isa(this.ParamKernel,'kernels.GaussKernel') && ~isa(this.ParamKernel,'kernels.NoKernel'))
                error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
            end    
            
            %% Initializations
            atd = model.Data.ApproxTrainData;
            % get validation data set vd
            vnum = round(size(atd,2)*this.ValidationPercent);
            sel = round(linspace(1,size(atd,2),vnum));
            vd = atd(:,sel);
            atd(:,sel) = [];
            
            fx = model.Data.ApproxfValues;
            vdfx = fx(:,sel);
            fx(:,sel) = [];
            
            % check if training data contains parameters
            if any(atd(1,:) ~= 0)
                hasparams = true;
                params = model.Data.getParams(atd(1,:));
                vdparams = model.Data.getParams(vd(1,:));
            else
                hasparams = false;
                params = [];
                vdparams = [];
            end
            nx = general.NNTracker;
            nt = general.NNTracker;
            if hasparams
                np = general.NNTracker;
            end
            
            %% Compute bounding boxes & center
            [BXmin, BXmax] = general.Utils.getBoundingBox(atd(4:end,:));
            Btmin = min(atd(3,:)); Btmax = max(atd(3,:));
            % Get bounding box diameters
            bxdia = norm(BXmax - BXmin);
            btdia = Btmax-Btmin;
            if hasparams
                [BPmin, BPmax] = general.Utils.getBoundingBox(...
                    model.Data.getParams(unique(atd(1,:))));
                thecenter = [BXmin+BXmax; Btmin+Btmax; BPmin + BPmax]/2;
                B = [atd(3:end,:); params];
                bpdia = norm(BPmax - BPmin);
            else
                thecenter = [BXmin+BXmax; Btmin+Btmax]/2;
                B = atd(3:end,:);
            end
            
            %% Select initial center x0
            % Strategy: Take the point that is closest to the bounding
            % box center!
            A = repmat(thecenter, 1, size(atd,2));
            [dummy, inIdx] = min(sum((A-B).^2,1));
            clear A B;
            
            initial = atd(:,inIdx);
            this.Centers.xi = initial(4:end);
            nx.addPoint(initial(4:end)); % Add points to nearest neighbor trackers (for gamma comp)
            this.Centers.ti = initial(3);
            nt.addPoint(initial(3));
            if hasparams
                this.Centers.mui = params(:,inIdx);
                np.addPoint(this.Centers.mui);
            else
                this.Centers.mui = [];
            end
            
            %% Set up initial expansion
            used = inIdx;
            this.Ma = fx(:,inIdx);
            this.off = [];
            
            %% Choose initial gammas
            dists = [dfun(minfac*bxdia, fac*bxdia); dfun(minfac*btdia, fac*btdia)];
            if hasparams
                dists = [dists; dfun(minfac*bpdia, fac*bpdia)];
            end
            minerr = Inf; gt = []; gp = [];
            for idx = 1:size(dists,2)
                gx = this.SystemKernel.setGammaForDistance(dists(1,idx),this.gameps);
                if ~isa(this.TimeKernel,'kernels.NoKernel')
                    gt = this.TimeKernel.setGammaForDistance(dists(2,idx),this.gameps);
                end
                if hasparams && ~isa(this.ParamKernel,'kernels.NoKernel')
                    gp = this.ParamKernel.setGammaForDistance(dists(3,idx),this.gameps);
                end
                
                fhat = this.evaluate(atd(4:end,:),atd(3,:),params);
                val = errfun(fx,fhat);
                if val < minerr
                    minerr = val;
                    bestgx = gx;
                    bestgt = gt;
                    bestgp = gp;
                end
            end
            %% Assign best values
            this.SystemKernel.Gamma = bestgx;
            if ~isa(this.TimeKernel,'kernels.NoKernel')
                this.TimeKernel.Gamma = bestgt;
            end
            if hasparams && ~isa(this.ParamKernel,'kernels.NoKernel')
                this.ParamKernel.Gamma = bestgp;
            end
            if KerMor.App.Verbose > 1
                fprintf('Initial gammas: SK:%e, TK:%e, PK:%e\n',bestgx,bestgt,bestgp);
            end
            
            %% Outer control loop
            cnt = 1;
            % Stopping condition preps
            this.effabs = this.MaxAbsErrFactor * max(abs(model.Data.ApproxfValues(:)));
            
            % Keep track of maximum errors
            this.MaxErrors = zeros(size(1,this.MaxExpansionSize));
            while true
                
                %% Determine maximum error over training data
                fhat = this.evaluate(atd(4:end,:),atd(3,:),params);
                [val, maxidx, errs] = errfun(fx,fhat);
                rel = val / (norm(model.Data.ApproxfValues(maxidx))+eps);
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
                new = atd(:,maxidx);
                this.Centers.xi(:,end+1) = new(4:end);
                this.Centers.ti(end+1) = new(3);
                if hasparams
                    this.Centers.mui(:,end+1) = params(:,maxidx);
                    np.addPoint(params(:,maxidx));
                end
                % Add points to nearest neighbor trackers (for gamma comp)
                nx.addPoint(new(4:end));
                nt.addPoint(new(3));
                
                %% Compute new approximation
                dists = [dfun(nx.getMinNN, fac*bxdia); dfun(nt.getMinNN, fac*btdia)];
                if hasparams
                    minnn = bpdia/this.NumGammas;
                    if ~isinf(np.getMinNN)
                        minnn = np.getMinNN;
                    end
                    dists = [dists; dfun(minnn, fac*bpdia)];%#ok
                end
                
                if KerMor.App.Verbose > 2
                    %dists%#ok
                end
                
                minerr = Inf; gt = []; gp = [];
                for gidx = 1:size(dists,2)
                    d = dists(:,gidx);
                    % Update Kernels Gamma values
                    gx = this.SystemKernel.setGammaForDistance(d(1),this.gameps);
                    %if KerMor.App.Verbose > 2
                        %fprintf('Kernels - Sys:%10f => gamma=%f',d(1),gx);
                        fprintf('xg:%.5e',gx);
                    %end
                    if ~isa(this.TimeKernel,'kernels.NoKernel')
                        gt = this.TimeKernel.setGammaForDistance(d(2),this.gameps);
                        %if KerMor.App.Verbose > 2
                            %fprintf(', Time:%10f => gamma=%10f',d(2),gt);
                            fprintf(', tg:%.5e',gx);
                        %end
                    end
                    if hasparams && ~isa(this.ParamKernel,'kernels.NoKernel')
                        gp = this.ParamKernel.setGammaForDistance(d(3),this.gameps);
                        %if KerMor.App.Verbose > 2
                            fprintf(', pg=%.5e',gp);
                        %end
                    end

                    %% Compute coefficients
                    warning('off','MATLAB:nearlySingularMatrix');
                    % Call coeffcomp preparation method and pass kernel matrix
                    K = this.getKernelMatrix;
                    this.CoeffComp.init(K);

                    % Call protected method
                    this.computeCoeffs(fx(:,used));
                    warning('on','MATLAB:nearlySingularMatrix');
                    
                    % get error on training data
                    fhat = this.evaluate(atd(4:end,:),atd(3,:),params);
                    val = errfun(fx, fhat);
                    impro = (val / minerr) * 100;
                    
                    % get error on validation set
                    fvali = this.evaluate(vd(4:end,:),vd(3,:),vdparams);
                    vdval = errfun(vdfx, fvali);
                    if val < minerr
                        minerr = val;
                        bestgx = gx;
                        bestgt = gt;
                        bestgp = gp;
                        bestMa = this.Ma;
                        bestoff = this.off;
                        %if KerMor.App.Verbose > 2
                            fprintf(' b: %.5e, %3.2f%%',val,impro);
                        %end
                    else
                        %if KerMor.App.Verbose > 2
                            %break;
                            fprintf(' w: %.5e, %3.2f%%',val,impro);
                        %end
                    end
                    
                    %if KerMor.App.Verbose > 2
                        fprintf(' ||Ma||:%.5e, vd-err:%.5e, prod:%.5e\n',sum(sqrt(sum(this.Ma.^2,1))),vdval,val*vdval);
                    %end
                
                end
                fprintf('\n');
                
                %% Assign best values
                this.SystemKernel.Gamma = bestgx;
                if ~isa(this.TimeKernel,'kernels.NoKernel')
                    this.TimeKernel.Gamma = bestgt;
                end
                if hasparams && ~isa(this.ParamKernel,'kernels.NoKernel')
                    this.ParamKernel.Gamma = bestgp;
                end
                this.Ma = bestMa;
                this.off = bestoff;
                
                %if KerMor.App.Verbose > 1
                    fprintf('-- It: %d ---- Minerr: %f ----- Best values: System:%f, Time:%f, Param:%f ----------\n',cnt,minerr,bestgx,bestgt,bestgp);
                %end
                
                cnt = cnt+1;
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
            
            function d = logsp(from, to)%# ok
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
                    if size(this.Centers.xi,2) > 1
                        plot(this.Centers.xi(:,1:end-1),'.','MarkerSize',2);
                    end
                    plot(this.Centers.xi(:,end),'r*','MarkerSize',3);
                    
                    % Plot params & time also
                    if KerMor.App.Verbose > 3
                        subplot(2,2,3); hold off;
                        plot([Btmin, Btmax],'black');
                        axis tight;
                        hold on;
                        if size(this.Centers.ti,2) > 1
                            plot(this.Centers.ti(:,1:end-1),'.','MarkerSize',2);
                        end
                        plot(this.Centers.ti(:,end),'r*','MarkerSize',3);
                        subplot(2,2,4); hold off;
                        plot([BPmin, BPmax],'black');
                        axis tight;
                        hold on;
                        if size(this.Centers.mui,2) > 1
                            plot(this.Centers.mui(:,1:end-1),'.','MarkerSize',2);
                        end
                        plot(this.Centers.mui(:,end),'r*','MarkerSize',3);
                    end
                    
                   pause;
            end
     
        end
                        
        function target = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            target = approx.AdaptiveCompWiseKernelApprox;
            
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
            
            % copy local props
            copy.MaxExpansionSize = this.MaxExpansionSize;
            copy.NumGammas = this.NumGammas;
            %copy.dfact = this.dfact;
            copy.gameps = this.gameps;
            copy.MaxRelErr = this.MaxRelErr;
            copy.MaxAbsErrFactor = this.MaxAbsErrFactor;
            copy.MaxErrors = this.MaxErrors;
        end       
        
        function set.MaxExpansionSize(this, value)
            if ~isposintscalar(value)
                error('Value must be a positive integer.');
            end
            this.MaxExpansionSize = value;
        end
        
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
        function bool = checkStop(this,cnt,rel,val)
            % Checks the stopping conditions 
            bool = false;
            if cnt == this.MaxExpansionSize
                disp('Stopping reason: Max expansion size reached.');
                bool = true;
            elseif rel < this.MaxRelErr
                fprintf('Stopping reason: relative error %.7e < %.7e\n',rel,this.MaxRelErr);
                bool = true;
            elseif val < this.effabs
                fprintf('Stopping reason: absolute error %.7e < %.7e\n',val,this.effabs);
                bool = true;
            end
        end
    end
end


