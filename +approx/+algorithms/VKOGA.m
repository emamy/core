classdef VKOGA < approx.algorithms.AAdaptiveBase
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
        UsefPGreedy = false;
        
        MaxAbsResidualErr = 1e-5;
        
        % Determines which error measure is to use to select a solution if the algorithm stops
        % for any reason but matching one of the MaxRelErr or MaxAbsResidualErr tolerances.
        %
        % Allowed values are 'abs' and 'rel'
        %
        % @type char @default 'abs'
        FailureErrorMeasure = 'abs';
        
%         Gain;
%         HerrDecay;
% 
%         % Lower bound for orthonormal remainders (kernel matrix conditioning improvement for
%         % larger bounds
%         PhiNormMin = sqrt(eps);
%         
%         % The original kernel expansion used to generate the approximation
%         % training data (DEBUG)
%         f;
%         
        VKOGABound;
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
%             copy.Gain = this.Gain;
%             copy.HerrDecay = this.HerrDecay;
            copy.MaxRelErrors = this.MaxRelErrors;
            copy.FailureErrorMeasure = this.FailureErrorMeasure;
%             copy.used = this.used;
%             copy.f = this.f;
            copy.VKOGABound = this.VKOGABound;
%             copy.PhiNormMin = this.PhiNormMin;
        end
    end
    
    methods(Access=protected, Sealed)
        function startAdaptiveExtension(this, kexp, atd)
            % Starts the adaptive extension of the VKOGA algorithm.
            
            % Flag for experimental mode
            exp_mode = 1 == 0;
            vb = KerMor.App.Verbose;
            
            ec = this.ExpConfig;
            nc = ec.getNumConfigurations;
            
            this.basis_norms = this.MaxErrors;
            this.StopFlags = zeros(nc,1);
            this.VKOGABound = this.basis_norms;

            xi = atd.xi.toMemoryMatrix;
            fxi = atd.fxi.toMemoryMatrix;
            fxinorm = this.ErrorFun(fxi);
            fxinorm(fxinorm == 0) = 1;
            N = size(xi,2);
            kexp.clear;
            
            minerr = Inf;
            minrelerr = Inf;
            bestNV = [];
            bestc = [];

            if exp_mode
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            
            % Kernel matrix diagonal
            Kdiag = ones(1,N);
            
            %% Run loop for all desired distances
            pi = ProcessIndicator('VKOGA approximation for %d kernel configurations',nc,false,nc);
            for cidx = 1:nc
                m = 1;
                
                % Set current hyperconfiguration
                ec.applyConfiguration(cidx, kexp);
                
                % Values of Newton basis
                NV = zeros(N,this.MaxExpansionSize);
                NV(:,1) = kexp.getKernelMatrixColumn(this.initialidx, xi);
                
                % Coefficients of Newton basis
                c = zeros(size(atd.fxi,1),this.MaxExpansionSize);
                c(:,1) = atd.fxi(:,this.initialidx)/sqrt(Kdiag(this.initialidx));
                
                % Sum_j N_j^2(x_i) (for denominator)
                sumNsq = (NV(:,1).^2)';
                
                fresidual = fxi;
                
                free = true(1,N);
                free(this.initialidx) = false;
                used = zeros(1,N);
                used(1) = this.initialidx;
                
                this.VKOGABound(cidx, m) = 1/max(Kdiag - sumNsq);
                               
                %% Main extension loop
                while true
                                        
                    %% Residual and error computation                    
                    fresidual = fresidual - c(:,m)*(NV(:,m))'; % Cumulative computation (eff.)
%                     fresidual1 = fxi - c(:,1:m)*(NV(:,1:m))'; % Complete computation
%                     figure(1);
%                     semilogy(abs(fresidual1 - fresidual));

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
                        div(div <= 0) = Inf;
                    else
                        div = 1;
                    end
                    sum_fresidual = sum(fresidual.^2,1) ./ div;
                    
                    [~, maxidx] = max(sum_fresidual);
                    tN = kexp.getKernelMatrixColumn(maxidx, xi) - NV(:,1:m)*NV(maxidx,1:m)';
                    
                    if exp_mode
                        pm.resetCount;
                        h1 = pm.nextPlot('fun','Function','x','f(x)');
                        plot(h1,xi(1,:),[fxi; c(:,1:m)*(NV(:,1:m))']');
                        h2 = pm.nextPlot('nfun','Newton Basis Function','x','N_i(x)');
                        plot(h2,xi(1,:),NV(:,1:m)); 
                        h3 = pm.nextPlot('err','Absolute error','x','|f(x)-f^m(x)|');
                        LogPlot.cleverPlot(h3,1:m,this.MaxErrors(:,1:m)); 
                        h4 = pm.nextPlot('MaxRelErrors','Relative error','x','|(f(x)-f^m(x))/f(x)|');
                        LogPlot.cleverPlot(h4,1:m,this.MaxRelErrors(:,1:m));
                        pm.done;
                    end
                    
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
                    free(maxidx) = false;
                    used(m) = maxidx;
                    
                    % Debug stuff
                    this.basis_norms(cidx,m) = tNnorm;
                    this.VKOGABound(cidx, m) = this.VKOGABound(cidx, m-1) + 1/max(Kdiag - sumNsq);
%                     if exp_mode
%                         h = pm.nextPlot('sumNSq','Power fun','x','P(x)');
%                         LogPlot.cleverPlot(h,xi(1,:),sumNsq); 
%                         pm.done;
%                     end
                end
                
                better = false;
                if stopflag == StopFlag.ABS_ERROR
                    bestm = m;
                    if this.MaxErrors(cidx,m) < minerr
                        better = true;
                        minerr = this.MaxErrors(cidx,m);
                    end
                elseif stopflag == StopFlag.REL_ERROR
                    bestm = m;
                    if this.MaxRelErrors(cidx,m) < minrelerr
                        better = true;
                        minrelerr = this.MaxRelErrors(cidx,m);
                    end
                elseif stopflag == StopFlag.NEGATIVE_POWFUN || stopflag == StopFlag.MAX_SIZE
                    if strcmp(this.FailureErrorMeasure,'rel')
                        [val, bestm] = min(this.MaxRelErrors(cidx,1:m));
                        if val < minrelerr
                            better = true;
                            minrelerr = val;
                        end
                    else
                        [val, bestm] = min(this.MaxErrors(cidx,1:m));
                        if val < minerr
                            better = true;
                            minerr = val;
                        end
                    end
                end
                this.ExpansionSizes(cidx) = bestm;
                if better
                    bestused = used(1:bestm);
                    bestNV = NV(:,1:bestm);
                    bestc = c(:,1:bestm);
                    bestcidx = cidx;
                end
                this.StopFlags(cidx) = stopflag;
                pi.step;
            end
            
            %% Set (& apply) best configuration
            this.Used = bestused;
            ec.setBestConfig(bestcidx, kexp);
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
                    bestcidx,og,length(bestused));
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
                if isfield(this,'UsefPGreedy')
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
            % Tests the VKOGA algorithm
            
            demos.VKOGA.demo_VKOGA_1D_nD(1,false);
            demos.VKOGA.demo_VKOGA_1D_nD(1,true);
            demos.VKOGA.demo_VKOGA_1D_nD(2,false);
            demos.VKOGA.demo_VKOGA_1D_nD(2,true);
            demos.VKOGA.demo_VKOGA_1D_nD(3,false,10);
            demos.VKOGA.demo_VKOGA_1D_nD(3,true);
            demos.VKOGA.demo_VKOGA_1D_nD(4,false,10);
            demos.VKOGA.demo_VKOGA_1D_nD(4,true);
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
            
            kexp = kernels.KernelExpansion;
            k = kernels.GaussKernel(.8);
            kexp.Kernel = k;
            atd = data.ApproxTrainData(xi,[],[]);
            atd.fxi = fxi;
            
            %% Algorithm setup
            alg = approx.algorithms.VKOGA;
            alg.MaxExpansionSize = 100;
            alg.UsefPGreedy = true;
            ec = kernels.config.ExpansionConfig;            
            gammas = linspace(.5,2,5);
            ec.StateConfig = kernels.config.RBFConfig('G',gammas);
            alg.ExpConfig = ec;

            alg.computeApproximation(kexp, atd);
            
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
            
            m = length(alg.Used);
            h = pm.nextPlot('err','Absolute error','x','|f(x)-f^m(x)|');
            ph = LogPlot.cleverPlot(h,1:m,alg.MaxErrors(:,1:m));
            set(ph(alg.ExpConfig.vBestConfigIndex),'LineWidth',2);
            
            h = pm.nextPlot('MaxRelErrors','Relative error','x','|(f(x)-f^m(x))/f(x)|');
            ph = LogPlot.cleverPlot(h,1:m,alg.MaxRelErrors(:,1:m));
            set(ph(alg.ExpConfig.vBestConfigIndex),'LineWidth',2);
            pm.done;
            
            d = struct;
            d.kexp = kexp;
            d.atd = atd;
            d.alg = alg;
            res = true;
        end
    end
end
