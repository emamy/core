classdef FixedCompWiseKernelApprox < approx.algorithms.BaseKernelApproxAlgorithm
% FixedCompWiseKernelApprox: Component-wise kernel approximation with fixed center set.
%
% @author Daniel Wirtz @date 2011-06-01
%
% @change{0,7,dw,2011-11-22} Moved FixedCompWiseKernelApprox.guessGammas to
% kernels.config.ExpansionConfig
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
        % The different kernel expansion configurations to try
        %
        % @propclass{critical} Without this setting this algorithm makes little sense.
        %
        % @type kernels.config.ExpansionConfig @default []
        ExpConfig = [];
        
        % The different coefficient computation algorithm configurations to try
        %
        % @propclass{critical} Without this setting this algorithm makes little sense.
        %
        % @type general.config.IClassConfig @default []
        CoeffConfig = [];
    end
    
    properties(SetAccess=private)
        BestConfigIndices;
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
            this.registerProps('ExpConfig');
        end
                        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.FixedCompWiseKernelApprox;
            
            copy = clone@approx.algorithms.BaseKernelApproxAlgorithm(this, copy);

            % copy local props
            copy.ExpConfig = this.ExpConfig;
        end
    end
    
    methods(Access=protected, Sealed)
        function detailedComputeApproximation(this, kexp, atd)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
             
            % Set centers
            xi = atd.xi.toMemoryMatrix;
            fxi = atd.fxi.toMemoryMatrix;
            kexp.Centers.xi = xi;
            kexp.Centers.ti = atd.ti;
            kexp.Centers.mui = atd.mui;
            
            errfun = @getLInftyErr;
            ec = this.ExpConfig;
            nc = ec.getNumConfigurations;
            cc = this.CoeffConfig;
            nco = cc.getNumConfigurations;
            
            total = max(nc,1)*max(nco,1);
            % Keep track of maximum errors
            this.MaxErrors = zeros(1,total);
            minerr = Inf;
            bestcidx = [];
            bestcoidx = [];
            bestMa = [];
            pi = tools.ProcessIndicator('Trying %d configurations (%d kernel, %d coeffcomp)',...
                total,false,total,nc,nco);
            cnt = 1;
            for kcidx = 1:nc
                if KerMor.App.Verbose > 2
                    fprintf('Applying expansion config %s\n',ec.getConfigurationString(kcidx));
                end
                ec.applyConfiguration(kcidx, kexp);
                
                % Call coeffcomp preparation method and pass kernel matrix
                K = data.MemoryKernelMatrix(kexp.getKernelMatrix);
                this.CoeffComp.init(K, kexp);

                for coidx = 1:nco
                    if KerMor.App.Verbose > 2
                        fprintf('Applying %s config %s\n',class(this.CoeffComp),cc.getConfigurationString(coidx));
                    end
                    cc.applyConfiguration(coidx, this.CoeffComp);
                    
                    % Call protected method
                    this.computeCoeffs(kexp, atd.fxi, []);
%                     this.computeCoeffs(kexp, atd.fxi, kexp.Ma);
                    
                    % Determine maximum error
                    fhat = kexp.evaluate(xi, atd.ti, atd.mui);
                    [val, maxidx, errs] = errfun(fxi,fhat);
                    rel = val / (norm(fxi(maxidx))+eps);
                    this.MaxErrors(cnt) = val;
                    
                    if val < minerr
                        minerr = val;
                        bestMa = kexp.Ma;
                        bestcidx = kcidx;
                        bestcoidx = coidx;
                        if KerMor.App.Verbose > 3
                            fprintf(' b: %.5e (rel %.5e)',val,rel);
                        end
                    else
                        if KerMor.App.Verbose > 3
                            fprintf(' w: %.5e  (rel %.5e)',val,rel);
                        end
                    end

                    if KerMor.App.Verbose > 3
                        fprintf(' ||Ma||:%.5e\n',sum(kexp.Ma_norms));
                    end
                    pi.step;
                    cnt = cnt+1;
                end
            end
            
            %% Assign best values
            this.BestConfigIndices = [bestcidx; bestcoidx];
            ec.applyConfiguration(bestcidx, kexp);
            kexp.Ma = bestMa;
            
            if KerMor.App.Verbose > 1
                figure;
                plot(this.MaxErrors,'r');
            end
            
            pi.stop;
        
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


