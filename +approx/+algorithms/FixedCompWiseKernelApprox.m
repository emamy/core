classdef FixedCompWiseKernelApprox < approx.algorithms.BaseKernelApproxAlgorithm
% FixedCompWiseKernelApprox: Component-wise kernel approximation with fixed center set.
%
% @author Daniel Wirtz @date 2011-06-01
%
% @change{0,5,dw,2011-11-02} 
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
% - Re-enabled the clone method
% - Using the new data.ApproxTrainData for guessGammas
%
% @new{0,5,dw,2011-10-14}
% - Added the method
% FixedCompWiseKernelApprox.guessGammas as a helper
% method to determine suitable Gaussian kernel configurations. The
% algorithm used is basically the same as in AdaptiveCompWiseKernelApprox.
% - This algorithm now also works with kernels.ParamTimeKernelExpansion's
% 
% @new{0,5,dw,2011-07-07} Moved the old approx.FixedCompWiseKernelApprox class to this class.
%
% @new{0,4,dw,2011-06-01} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
%
% @todo think of new structure on how to combine this with the getDists in BaseAdaptiveCWKA (or
% extract convenience methods further, general concept of "KernelConfig")
    
    properties(SetObservable)    
        % @propclass{experimental} 
        Gammas;
    end
    
    properties(Transient, SetAccess=private)
        % Contains the maximum errors for each iteration/center extension step performed by the last
        % run of this algorithm.
        MaxErrors = [];
    end
        
    methods    
        function this = FixedCompWiseKernelApprox
            this = this@approx.algorithms.BaseKernelApproxAlgorithm;
            
            % Register default property changed listeners
            this.registerProps('Gammas');
        end
                        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.FixedCompWiseKernelApprox;
            
            copy = clone@approx.algorithms.BaseKernelApproxAlgorithm(this, copy);

            % copy local props
            copy.Gammas = this.Gammas;
        end

        function guessGammas(this, atd, ng, mf, Mf)
            % Guesses 'ng' `\gamma` values to use during approximation.
            %
            % The computation is using the same method as in BaseAdaptiveCWKA.getDists
            %
            % Parameters:
            % atd: The approximation training data @type data.ApproxTrainData
            % ng: The number of `\gamma` values to guess @type integer
            % mf: The factor for minimum sample data distance @type double @default .5
            % Mf: The factor for the bounding box diameter @type double @default 2
            %
            
            if ~isa(atd,'data.ApproxTrainData')
                error('First argument has to be a data.ApproxTrainData instance');
            end
            gameps = .1;
            if nargin < 5
                Mf = 2;
                if nargin < 4
                    mf = .5;
                end
            end
            dfun = @logsp;
            
            dists = zeros(3,ng);
            dists(1,:) = dfun(mf*atd.xiDia, Mf*atd.xiDia);
            if atd.hasTime
                dists(2,:) = dfun(mf*atd.tiDia, Mf*atd.tiDia);
            end
            if atd.hasParams
                dists(3,:) = dfun(mf*atd.muiDia, Mf*atd.muiDia);
            end
            k = kernels.GaussKernel;
            this.Gammas = zeros(size(dists));
            for i=1:ng
               this.Gammas(1,i) = k.setGammaForDistance(dists(1,i),gameps);
               this.Gammas(2,i) = k.setGammaForDistance(dists(2,i),gameps);
               this.Gammas(3,i) = k.setGammaForDistance(dists(3,i),gameps);
            end
            
            function linsp(from, to)%#ok
                d = linspace(from,to,ng);
            end
            
            function d = logsp(from, to)
                d = logspace(log10(from),log10(to),ng);
            end
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedComputeApproximation(this, kexp, atd)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            pte = isa(kexp,'kernels.ParamTimeKernelExpansion');
            if pte
                tkg = isa(kexp.TimeKernel,'kernels.GaussKernel');
                if ~isa(kexp.Kernel,'kernels.GaussKernel') || ...
                        (pte && ((~tkg && ~isa(kexp.TimeKernel,'kernels.NoKernel')) || ...
                        ~isa(kexp.ParamKernel,'kernels.GaussKernel')))
                    error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
                end
            end
            
            % Set AKernelCoreFun centers
            kexp.Centers.xi = atd.xi;
            kexp.Centers.ti = atd.ti;
            kexp.Centers.mui = atd.mui;
            
            errfun = @getLInftyErr;
                        
            % Keep track of maximum errors
            this.MaxErrors = zeros(1,length(this.Gammas));
            minerr = Inf;
            bestg = []; 
            bestMa = [];
            for gidx = 1:size(this.Gammas,2)
                
                g = this.Gammas(:,gidx);
                kexp.Kernel.Gamma = g(1);
                if pte
                    if tkg
                        kexp.TimeKernel.Gamma = g(2);
                    end
                    kexp.ParamKernel.Gamma = g(3);
                end
                
                if KerMor.App.Verbose > 2
                    fprintf('xg: %.5e',g(1));
                    if pte
                        fprintf(', tg: %.5e, pg:%.5e',g(1),g(2));
                    end
                    fprintf('\n');
                end
                
                %% Compute coefficients
                %warning('off','MATLAB:nearlySingularMatrix');
                % Call coeffcomp preparation method and pass kernel matrix
                K = data.MemoryKernelMatrix(kexp.getKernelMatrix);
                this.CoeffComp.init(K, kexp);

                % Call protected method
                ex = [];
                try
                    this.computeCoeffs(kexp, atd.fxi, []);
                catch ME
                    if gidx > 1 && strcmp(ME.identifier,'KerMor:coeffcomp:failed')    
                        ex = ME;
                        break;
                    else
                        rethrow(ME);
                    end
                end
                
                %% Determine maximum error over training data
                fhat = kexp.evaluate(atd.xi, atd.ti, atd.mui);
                [val, maxidx, errs] = errfun(atd.fxi,fhat);
                rel = val / (norm(atd.fxi(maxidx))+eps);
                this.MaxErrors(gidx) = val;
                
                if val < minerr
                    minerr = val;
                    bestMa = kexp.Ma;
                    bestg = g;
                    if KerMor.App.Verbose > 2
                        fprintf(' b: %.5e (rel %.5e)',val,rel);
                    end
                else
                    if KerMor.App.Verbose > 2
                        fprintf(' w: %.5e  (rel %.5e)',val,rel);
                    end
                end
                    
                if KerMor.App.Verbose > 2
                    fprintf(' ||Ma||:%.5e\n',sum(kexp.Ma_norms));
                end
            end
            
                
            %% Assign best values
            kexp.Kernel.Gamma = bestg(1);
            if pte
                if isa(kexp.TimeKernel,'kernels.GaussKernel')
                    kexp.TimeKernel.Gamma = bestg(2);
                end
                kexp.ParamKernel.Gamma = bestg(3);
            end
            kexp.Ma = bestMa;

            if ~isempty(ex)
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
        end
    end
end


