classdef VKOGA < approx.algorithms.BaseAdaptiveCWKA
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
        err;
        relerr;
        expsizes;
        
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
%         VKOGABound;
    end
    
    properties(SetAccess=private)
        bestNewtonBasisValuesOnATD;
        
        % debug props
        basis_norms;
    end
    
    methods
        function this = VKOGA
            this = this@approx.algorithms.BaseAdaptiveCWKA;
        end

        function copy = clone(this)
            % Clones the instance.
            copy = approx.algorithms.VectorialKernelOMP;
            copy = clone@approx.algorithms.BaseAdaptiveCWKA(this, copy);
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
        
        function plotErrors(this, pm)
            if nargin < 2
                pm = PlotManager(false,1,2);
                pm.LeaveOpen = true;
            end
            
            h = pm.nextPlot('abs','Absolute errors','step','value');
            ph = semilogy(h,1:size(this.err,2),this.err');
            set(ph(this.ExpConfig.vBestConfigIndex),'LineWidth',2);
            h = pm.nextPlot('rel','Relative errors','step','value');
            ph = semilogy(h,1:size(this.relerr,2),this.relerr');
            set(ph(this.ExpConfig.vBestConfigIndex),'LineWidth',2);
            
            if nargin < 2
                pm.done;
            end
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
            
            this.err = zeros(nc,this.MaxExpansionSize);
            this.relerr = this.err;
            this.basis_norms = this.err;
            this.expsizes = zeros(nc,1);

            xi = atd.xi.toMemoryMatrix;
            fxi = atd.fxi.toMemoryMatrix;
            fxinorm = this.ErrorFun(fxi);
            fxinorm(fxinorm == 0) = 1;
            N = size(xi,2);
            kexp.clear;
            
            minerr = Inf;
            bestNV = [];
            bestc = [];

            if exp_mode
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            
            %% Run loop for all desired distances
            pi = tools.ProcessIndicator('VKOGA approximation for %d gamma values',nc,false,nc);
            for cidx = 1:nc
                m = 1;
                
                % Set current hyperconfiguration
                ec.applyConfiguration(cidx, kexp);
                
                % Kernel matrix diagonal
                Kdiag = ones(1,N);
                
                % Values of Newton basis
                NV = zeros(N,this.MaxExpansionSize);
                NV(:,1) = kexp.getKernelMatrixColumn(this.initialidx, xi);
                
                % Coefficients of Newton basis
                c = zeros(size(atd.fxi,1),this.MaxExpansionSize);
                c(:,1) = atd.fxi(:,this.initialidx)/sqrt(Kdiag(this.initialidx));
                
                % Sum_j N_j^2(x_i) (for denominator)
                sumNsq = (NV(:,1).^2)';
                
                cum_fresidual = fxi;
                
                free = true(1,N);
                free(this.initialidx) = false;
                used = zeros(1,N);
                used(1) = this.initialidx;
                               
                %% Main extension loop
                while true
                                        
                    %% Residual and error computation
                    fresidual = fxi - c(:,1:m)*(NV(:,1:m))';
                    cum_fresidual = cum_fresidual - c(:,m)*(NV(:,m))';
                    e = this.ErrorFun(cum_fresidual);
                    this.err(cidx,m) = max(e);
                    this.relerr(cidx,m) = max(e./fxinorm);
                    
                    %% Check stopping conditions
                    if m == this.MaxExpansionSize
                        if vb > 1
                            fprintf('VKOGA stopping criteria holds: Max expansion size %d reached.\nResidual error %.7e > %.7e, Max relative error %.7e > %.7e\n',...
                                m,this.err(cidx,m-1),this.MaxAbsResidualErr,this.relerr(cidx,m-1),this.MaxRelErr);
                        end
                        break;
                    elseif this.relerr(cidx,m) < this.MaxRelErr
                        if vb > 1
                            fprintf('VKOGA stopping criteria holds: Relative error %.7e < %.7e\n',rel,this.MaxRelErr);
                        end
                        break;
                    elseif this.err(cidx,m) < this.MaxAbsResidualErr
                        if vb > 1
                            fprintf('VKOGA stopping criteria holds: Residual error %.7e < %.7e\n',this.err(cidx,m),this.MaxAbsResidualErr);
                        end
                        break;
                    %% Emergency break: some 
                    elseif any(sumNsq(free) > 1)
                        if vb > 1
                            fprintf('VKOGA roundoff emergency stop at iteration %d: Negative power function values detected.\nResidual error %.7e > %.7e, Max relative error %.7e > %.7e\n',...
                                m,this.err(cidx,m-1),this.MaxAbsResidualErr,this.relerr(cidx,m-1),this.MaxRelErr);
                        end
                        break;
                    end
                    
                    %% Next basis point selection
                    if this.UsefPGreedy
                        % % Cap too small norms!
                        % phinormsq = this.PhiNormMin^2;
                        div = Kdiag - sumNsq;
                        div(div <= 0) = 1;
                    else
                        div = 1;
                    end
                    sum_fresidual = sum(fresidual.^2,1) ./ div;
                    
                    [~, maxidx] = max(sum_fresidual);
                    
                    if exp_mode
                        pm.resetCount;
                        h1 = pm.nextPlot('fun','Function','x','f(x)');
                        plot(h1,xi,[fxi; c(:,1:m)*(NV(:,1:m))']');
                        h2 = pm.nextPlot('nfun','Newton Basis Function','x','N_i(x)');
                        plot(h2,xi,NV(:,1:m)); 
                        h3 = pm.nextPlot('err','Absolute error','x','|f(x)-f^m(x)|');
                        tools.LogPlot.cleverPlot(h3,1:m,this.err(:,1:m)); 
                        h4 = pm.nextPlot('relerr','Relative error','x','|(f(x)-f^m(x))/f(x)|');
                        tools.LogPlot.cleverPlot(h4,1:m,this.relerr(:,1:m)); 
                        pm.done;
                    end
                    
                    tN = kexp.getKernelMatrixColumn(maxidx, xi) - NV(:,1:m)*NV(maxidx,1:m)';
                    m = m+1;
                    tNnorm = sqrt(Kdiag(maxidx)-sumNsq(maxidx));
                    NV(:,m) = tN/tNnorm;
                    c(:,m) = fresidual(:,maxidx)./tNnorm;
                    sumNsq = sumNsq + (NV(:,m).^2)';
                    this.basis_norms(cidx,m) = tNnorm;
                    free(maxidx) = false;
                    used(m) = maxidx;
                end
                
                [val, bestm] = min(this.err(cidx,1:m));
                this.expsizes(cidx) = bestm;
                if val < minerr
                    bestused = used(1:bestm);
                    bestNV = NV(:,1:bestm);
                    bestc = c(:,1:bestm);
                    minerr = val;
                    bestcidx = cidx;
                end
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
        
        function res = test_VKOGA2D1D
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
            ph = tools.LogPlot.cleverPlot(h,1:m,alg.err(:,1:m));
            set(ph(alg.ExpConfig.vBestConfigIndex),'LineWidth',2);
            
            h = pm.nextPlot('relerr','Relative error','x','|(f(x)-f^m(x))/f(x)|');
            ph = tools.LogPlot.cleverPlot(h,1:m,alg.relerr(:,1:m));
            set(ph(alg.ExpConfig.vBestConfigIndex),'LineWidth',2);
            pm.done;
            
            res.kexp = kexp;
            res.atd = atd;
            res.alg = alg;
        end
    end
end
