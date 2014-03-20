classdef VKOGA < approx.algorithms.AAdaptiveBase
% VKOGA: Vectorial kernel orthogonal greedy algorithm
%
% The VKOGA algorithm in introduced in @cite WH13 and adopts principles of
% greedy algorithms @cite LT05 @cite T08
%
% @author Daniel Wirtz @date 2012-02-09
%
% @new{0,7,dw,2014-01-24} Now can finally also use non-normalized kernels
% with the VKOGA algorithm.
%
% @change{0,7,dw,2012-11-26} Renamed to VKOGA and starting to build in
% IClassConfig interfaces
%
% @new{0,6,dw,2012-02-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % Flag on whether to use the f/P-Greedy variant instead of
        % f-Greedy. @default false @type logical
        UsefPGreedy = false;
        
        % Stopping criteria for the VKOGA algorithm.
        % The execution for on expansion configuration is stopped if the
        % maximum absolute pointwise L2-error is below this value.
        % @default 1e-5 @type double
        MaxAbsResidualErr = 1e-5;
        
        % Determines which error measure is to use to select the "best"
        % solution w.r.t. to the validation data set (the training data set
        % is used if not given)
        %
        % Allowed values are 'abs' and 'rel'
        %
        % @type char @default 'abs'
        ValidationErrorMeasure = 'abs';
        
        VKOGABound;
        
        AllExpansions;
    end
    
    properties(SetAccess=private)
        bestNewtonBasisValuesOnATD;
        
        % debug props
        basis_norms;
    end
    
    methods
        function this = VKOGA
            this = this@approx.algorithms.AAdaptiveBase;
        end

        function copy = clone(this)
            % Clones the instance.
            copy = approx.algorithms.VKOGA;
            copy = clone@approx.algorithms.AAdaptiveBase(this, copy);
            copy.UsefPGreedy = this.UsefPGreedy;
            copy.MaxAbsResidualErr = this.MaxAbsResidualErr;
            copy.MaxRelErrors = this.MaxRelErrors;
            copy.ValidationErrorMeasure = this.ValidationErrorMeasure;
            copy.VKOGABound = this.VKOGABound;
            copy.basis_norms = this.basis_norms;
            copy.AllExpansions = this.AllExpansions;
        end
    end
    
    methods(Access=protected, Sealed)
        function kexp = startAdaptiveExtension(this, atd, avd)
            % Starts the adaptive extension of the VKOGA algorithm.
            vb = KerMor.App.Verbose;
            
            ec = this.ExpConfig;
            nc = ec.getNumConfigurations;
            
            this.basis_norms = this.MaxErrors;
            this.StopFlags = zeros(nc,1);
            this.VKOGABound = this.basis_norms;
            this.BestExpConfig = 0;
            this.AllExpansions = eval(sprintf('%s.empty',class(ec.Prototype)));

            xi = atd.xi.toMemoryMatrix;
            fxi = atd.fxi.toMemoryMatrix;
            fxinorm = this.ErrorFun(fxi);
            fxinorm(fxinorm == 0) = 1;
            if ~isempty(avd)
                vxi = avd.xi.toMemoryMatrix;
                vfxi = avd.fxi.toMemoryMatrix;
                vfxinorm = this.ErrorFun(vfxi);
                vfxinorm(vfxinorm == 0) = 1;
            end
            N = size(xi,2);
            
            minerr = Inf;
            bestNV = [];
            bestc = [];
            bestused = [];
            
            %% Run loop for all desired distances
            pi = ProcessIndicator('VKOGA approximation for %d kernel configurations',nc,false,nc);
            for cidx = 1:nc
                m = 1;
                
                % Set current hyperconfiguration
                kexp = ec.configureInstance(cidx);
                
                % Compute Kernel matrix diagonal
                if isa(kexp.Kernel,'kernels.ARBFKernel')
                    Kdiag = kexp.Kernel.evaluate(xi(:,1),xi(:,1))*ones(1,N);
                else
                    Kdiag = zeros(1,N);
                    for k=1:N
                        Kdiag(k) = kexp.Kernel.evaluate(xi(k,1),xi(k,1));
                    end
                end
                
                % Values of Newton basis
                NV = zeros(N,this.MaxExpansionSize);
                NV(:,1) = kexp.getKernelMatrixColumn(this.initialidx, xi);
                
                % Coefficients of Newton basis
                c = zeros(size(atd.fxi,1),this.MaxExpansionSize);
                c(:,1) = atd.fxi(:,this.initialidx)/sqrt(Kdiag(this.initialidx));
                
                % Sum_j N_j^2(x_i) (for denominator)
                sumNsq = (NV(:,1).^2)';
                
                fresidual = fxi;
                
                used = zeros(1,N);
                used(1) = this.initialidx;
                
                this.VKOGABound(cidx, m) = 1/max(Kdiag - sumNsq);
                               
                %% Main extension loop
                while true
                                        
                    %% Residual and error computation          
                    % Cumulative computation (more effective)
                    fresidual = fresidual - c(:,m)*(NV(:,m))'; 
                    % Complete computation
%                     fresidual1 = fxi - c(:,1:m)*(NV(:,1:m))'; 

                    e = this.ErrorFun(fresidual);
                    this.MaxErrors(cidx,m) = max(e);
                    this.MaxRelErrors(cidx,m) = max(e./fxinorm);
                    
                    %% Check stopping conditions
                    if m == this.MaxExpansionSize
                        if vb > 1
                            fprintf('VKOGA stopping criteria holds: Max expansion size %d reached.\nResidual error %.7e > %.7e, Max relative error %.7e > %.7e\n',...
                                m,this.MaxErrors(cidx,m-1),this.MaxAbsResidualErr,this.MaxRelErrors(cidx,m-1),this.MaxRelErr);
                        end
                        stopflag = StopFlag.MAX_SIZE;
                        break;
                    elseif this.MaxRelErrors(cidx,m) < this.MaxRelErr
                        if vb > 1
                            fprintf('VKOGA stopping criteria holds: Relative error %.7e < %.7e\n',this.MaxRelErrors(cidx,m),this.MaxRelErr);
                        end
                        stopflag = StopFlag.REL_ERROR;
                        break;
                    elseif this.MaxErrors(cidx,m) < this.MaxAbsResidualErr
                        if vb > 1
                            fprintf('VKOGA stopping criteria holds: Residual error %.7e < %.7e\n',this.MaxErrors(cidx,m),this.MaxAbsResidualErr);
                        end
                        stopflag = StopFlag.ABS_ERROR;
                        break;
                    end
                    
                    %% Next basis point selection
                    if this.UsefPGreedy
                        % Cap too small norms!
                        div = Kdiag - sumNsq;
                        div(used(1:m)) = Inf;
                        %div(div <= 0) = Inf;
                    else
                        div = 1;
                    end
                    sum_fresidual = sum(fresidual.^2,1) ./ div;
                    
                    [~, maxidx] = max(sum_fresidual);
                    tN = kexp.getKernelMatrixColumn(maxidx, xi) - NV(:,1:m)*NV(maxidx,1:m)';
                    
                    %% Emergency break:
                    % An entry of the power function for which a value greater than one has
                    % been summed up would be chosen, leading to imaginary norms.
                    if sumNsq(maxidx) > 1
                        if vb > 1
                            fprintf('VKOGA emergency stop at iteration %d: Power function value > 1 detected.\nResidual error %.7e > %.7e, Max relative error %.7e > %.7e\n',...
                                m,this.MaxErrors(cidx,m-1),this.MaxAbsResidualErr,this.MaxRelErrors(cidx,m-1),this.MaxRelErr);
                        end
                        stopflag = StopFlag.NEGATIVE_POWFUN;
                        break;
                    end
                    
                    m = m+1;
                    tNnorm = sqrt(Kdiag(maxidx)-sumNsq(maxidx));
                    NV(:,m) = tN/tNnorm;
                    c(:,m) = fresidual(:,maxidx)./tNnorm;
                    sumNsq = sumNsq + (NV(:,m).^2)';
                    used(m) = maxidx;
                    
                    this.basis_norms(cidx,m) = tNnorm;
                    this.VKOGABound(cidx, m) = this.VKOGABound(cidx, m-1) + 1/max(Kdiag - sumNsq);
                end
                
                this.ExpansionSizes(cidx) = m;
                this.MaxErrors(cidx,m+1:end) = Inf;
                this.MaxRelErrors(cidx,m+1:end) = Inf;
                
                kexp.setCentersFromATD(atd, used(1:m));
                kexp.Ma = c(:,1:m);
                kexp.Base = NV(used(1:m),1:m);
                this.AllExpansions(cidx) = kexp;
                
                % If set, compute error on validation set
                if ~isempty(avd)
                    e = this.ErrorFun(vfxi - kexp.evaluate(vxi));
                    this.ValidationErrors(cidx,1) = max(e);
                    this.ValidationRelErrors(cidx,1) = max(e./vfxinorm);
                    this.ValidationErrors(cidx,2) = mean(e);
                    this.ValidationRelErrors(cidx,2) = mean(e./vfxinorm);
                else
                    % Otherwise just copy the training set error
                    this.ValidationErrors(cidx,1) = this.MaxErrors(cidx,m);
                    this.ValidationRelErrors(cidx,1) = this.MaxRelErrors(cidx,m);
                    this.ValidationErrors(cidx,2) = mean(e);
                    this.ValidationRelErrors(cidx,2) = mean(e./fxinorm);
                end
                
                if strcmp(this.ValidationErrorMeasure,'rel')
                    val = this.ValidationRelErrors(cidx,2);
                else
                    val = this.ValidationErrors(cidx,2);
                end
                if val < minerr
                    minerr = val;
                    bestused = used(1:m);
                    bestNV = NV(:,1:m);
                    bestc = c(:,1:m);
                    this.BestExpConfig = cidx;
                end
                this.StopFlags(cidx) = stopflag;
                pi.step;
            end
            
            %% Set (& apply) best configuration
            if isempty(bestused)
                warning();
            end
            this.Used = bestused;
            kexp = ec.configureInstance(this.BestExpConfig);
            kexp.setCentersFromATD(atd, bestused);
            
            % Compute Ma
            this.bestNewtonBasisValuesOnATD = bestNV;
            kexp.Ma = bestc;
            % Extract partial cholesky matrix
            kexp.Base = bestNV(bestused,1:length(bestused));
            
            % Some cleanup
            maxsize = max(this.ExpansionSizes);
            this.MaxErrors = this.MaxErrors(:,1:maxsize);
            this.MaxRelErrors = this.MaxRelErrors(:,1:maxsize);
            
            % Debug stuff
            this.VKOGABound = size(fxi,1) ./ (1 + this.VKOGABound/size(fxi,1));
            
            if vb > 1
                og = 'off';
                if this.UsefPGreedy
                    og = 'on';
                end
                fprintf('VKOMP best kernel config index:%d for VKOGA with f/P-Greedy=%s, exp-size:%d\n',...
                    this.BestExpConfig,og,length(bestused));
            end
            
            pi.stop;
        end
    end
    
    methods(Static,Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'approx.algorithms.VKOGA')
                a = approx.algorithms.VKOGA;
                if isfield(this,'UsefPGreedy') && ~isempty(this.UsefPGreedy)
                    a.UsefPGreedy = this.UsefPGreedy;
                elseif isfield(this,'UseOGA')
                    a.UsefPGreedy = this.UseOGA;
                end
                if isfield(this,'MaxAbsResidualErr')
                    a.MaxAbsResidualErr = this.MaxAbsResidualErr;
                end
                if isfield(this,'bestNewtonBasisValuesOnATD')
                    a.bestNewtonBasisValuesOnATD = this.bestNewtonBasisValuesOnATD;
                end
                if isfield(this,'basis_norms')
                    a.basis_norms = this.basis_norms;
                end
                if isfield(this,'stopFlags')
                    a.StopFlags = this.stopFlags;
                end
                this = loadobj@approx.algorithms.AAdaptiveBase(a, this);
            end
            if isempty(this.StopFlags)
                this.StopFlags = zeros(size(this.MaxErrors));
            end
        end
    end
    
    methods(Static)
        
        function res = test_VKOGA1DnD
            % Tests the VKOGA algorithm by running the demos
            
            demos.VKOGA.VKOGA_1D_nD(1,false);
            demos.VKOGA.VKOGA_1D_nD(1,true);
            demos.VKOGA.VKOGA_1D_nD(2,false);
            demos.VKOGA.VKOGA_1D_nD(2,true);
            demos.VKOGA.VKOGA_1D_nD(3,false,10);
            demos.VKOGA.VKOGA_1D_nD(3,true);
            demos.VKOGA.VKOGA_1D_nD(4,false,10);
            demos.VKOGA.VKOGA_1D_nD(4,true);
            res = true;
        end
        
        function [res, d] = test_VKOGA2D1D
            % Tests the VKOGA algorithm
            
            %% Data setup
            [f, r] = testing.TestFunctions.F8;
            step = .02;
            x1 = r(1,1):step:r(1,2);
            x2 = r(2,1):step:r(2,2);
            [X,Y] = meshgrid(x1,x2);
            xi = [X(:)'; Y(:)'];
            fxi = f(xi);
            
            atd = data.ApproxTrainData(xi,[],[]);
            atd.fxi = fxi;
            
            %% Algorithm setup
            alg = approx.algorithms.VKOGA;
            alg.MaxExpansionSize = 100;
            alg.UsefPGreedy = true;
            ec = kernels.config.ExpansionConfig;
            ec.Prototype.Kernel = kernels.Wendland;
            c = Utils.createCombinations(linspace(.5,2,5),[2 3]);
            ec.StateConfig = kernels.config.WendlandConfig(...
                'G',c(1,:),'S',c(2,:),'Dim',1);
            alg.ExpConfig = ec;

            kexp = alg.computeApproximation(atd);
            
            %% Plot approximated function
            
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            
            h = pm.nextPlot('fun','Function','x','f(x)');
            Z = reshape(fxi,length(x1),[]);
            surf(h,X,Y,Z);
            
            h = pm.nextPlot('fun','Approximation','x','kexp(x)');
            aZ = reshape(kexp.evaluate(xi),length(x1),[]);
            surf(h,X,Y,aZ);
            
            h = pm.nextPlot('fun','Error','x','f(x)-kexp(x)');
            surf(h,X,Y,abs(Z-aZ));
            
            h = pm.nextPlot('nfun','Last Newton Basis Function','x','N(x)');
            Z = reshape(alg.bestNewtonBasisValuesOnATD(:,end),length(x1),[]);
            surf(h,X,Y,Z);
            
            alg.plotSummary(pm, 'test_VKOGA');
            pm.done;
            
            d = struct;
            d.kexp = kexp;
            d.atd = atd;
            d.alg = alg;
            res = true;
        end
    end
end
