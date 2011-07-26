classdef FixedCompWiseKernelApprox < approx.algorithms.BaseKernelApproxAlgorithm
% FixedCompWiseKernelApprox: Component-wise kernel approximation with fixed center set.
%
% @author Daniel Wirtz @date 2011-06-01
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
%             copy.Gammas = this.Gammas;
%             
%             copy.MaxErrors = this.MaxErrors;
%         end
    end
    
    methods(Access=protected, Sealed)
        function detailedComputeApproximation(this, kexp, xi, ti, mui, fxi)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            if ~isa(kexp.Kernel,'kernels.GaussKernel') || ...
                    (~isa(kexp.TimeKernel,'kernels.GaussKernel') && ~isa(kexp.TimeKernel,'kernels.NoKernel')) || ...
                    (~isa(kexp.ParamKernel,'kernels.GaussKernel') && ~isa(kexp.ParamKernel,'kernels.NoKernel'))
                error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
            end
            
            % Set AKernelCoreFun centers
            kexp.Centers.xi = xi;
            kexp.Centers.ti = ti;
            kexp.Centers.mui = mui;
            
            errfun = @getLInftyErr;
                        
            % Keep track of maximum errors
            this.MaxErrors = zeros(1,length(this.Gammas));
            minerr = Inf; bestg = 0; bestMa = [];
            for gidx = 1:length(this.Gammas)
                
                g = this.Gammas(gidx);
                kexp.Kernel.Gamma = g;
                
                if KerMor.App.Verbose > 2
                    fprintf('xg:%.5e',g);
                end
                
                %% Compute coefficients
                %warning('off','MATLAB:nearlySingularMatrix');
                % Call coeffcomp preparation method and pass kernel matrix
                K = kexp.getKernelMatrix;
                this.CoeffComp.init(K);

                % Call protected method
                ex = [];
                try
                    this.computeCoeffs(kexp, fxi);
                catch ME
                    if gidx > 1 && strcmp(ME.identifier,'KerMor:coeffcomp:failed')    
                        ex = ME;
                        break;
                    else
                        rethrow(ME);
                    end
                end
                
                %% Determine maximum error over training data
                fhat = kexp.evaluate(xi,ti,mui);
                [val, maxidx, errs] = errfun(fxi,fhat);
                rel = val / (norm(fxi(maxidx))+eps);
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
            kexp.Kernel.Gamma = bestg;
%             if ~isa(kexp.TimeKernel,'kernels.NoKernel')
%                 kexp.TimeKernel.Gamma = bestgt;
%             end
%             if hasparams && ~isa(kexp.ParamKernel,'kernels.NoKernel')
%                 kexp.ParamKernel.Gamma = bestgp;
%             end
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


