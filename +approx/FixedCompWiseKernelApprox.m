classdef FixedCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
% FixedCompWiseKernelApprox: Component-wise kernel approximation with fixed center set.
%
% @author Daniel Wirtz @date 2011-06-01
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
            this = this@approx.BaseCompWiseKernelApprox;
            
            % Register default property changed listeners
            this.registerProps('Gammas');
        end
                        
        function target = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            target = approx.AdaptiveCompWiseKernelApprox;
            
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
            
            %this.cloneLocalProps(target,mfilename('class'));
            % copy local props
            copy.Gammas = this.Gammas;
            
            copy.MaxErrors = this.MaxErrors;
        end
    end
    
    methods(Access=protected, Sealed)
        function computeCompwiseApprox(this, xi, ti, mui, fxi)
            % Performs adaptive approximation generation.
            %
            % @docupdate
            % @todo Think about suitable stopping condition (relative error
            % change?)
            
            %% Checks
            % This algorithm so far works only with Gaussian kernels
            if ~isa(this.SystemKernel,'kernels.GaussKernel') || ...
                    (~isa(this.TimeKernel,'kernels.GaussKernel') && ~isa(this.TimeKernel,'kernels.NoKernel')) || ...
                    (~isa(this.ParamKernel,'kernels.GaussKernel') && ~isa(this.ParamKernel,'kernels.NoKernel'))
                error('Any kernels used have to be Gaussian kernels for this approximation algorithm so far');
            end
            
            % Set AKernelCoreFun centers
            this.Centers.xi = xi;
            this.Centers.ti = ti;
            this.Centers.mui = mui;
            
            errfun = @getLInftyErr;
                        
            % Keep track of maximum errors
            this.MaxErrors = zeros(1,length(this.Gammas));
            minerr = Inf; bestg = 0; bestMa = [];
            for gidx = 1:length(this.Gammas)
                
                g = this.Gammas(gidx);
                this.SystemKernel.Gamma = g;
                
                if KerMor.App.Verbose > 2
                    fprintf('xg:%.5e',g);
                end
                
                %% Compute coefficients
                %warning('off','MATLAB:nearlySingularMatrix');
                % Call coeffcomp preparation method and pass kernel matrix
                K = this.getKernelMatrix;
                this.CoeffComp.init(K);

                % Call protected method
                ex = [];
                try
                    this.computeCoeffs(fxi);
                catch ME
                    if gidx > 1 && strcmp(ME.identifier,'KerMor:coeffcomp:failed')    
                        ex = ME;
                        break;
                    else
                        rethrow(ME);
                    end
                end
                
                %% Determine maximum error over training data
                fhat = this.evaluate(xi,ti,mui);
                [val, maxidx, errs] = errfun(fxi,fhat);
                rel = val / (norm(fxi(maxidx))+eps);
                this.MaxErrors(gidx) = val;
                
                if val < minerr
                    minerr = val;
                    bestMa = this.Ma;
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
                    fprintf(' ||Ma||:%.5e\n',sum(this.Ma_norms));
                end
            end
            
                
            %% Assign best values
            this.SystemKernel.Gamma = bestg;
%             if ~isa(this.TimeKernel,'kernels.NoKernel')
%                 this.TimeKernel.Gamma = bestgt;
%             end
%             if hasparams && ~isa(this.ParamKernel,'kernels.NoKernel')
%                 this.ParamKernel.Gamma = bestgp;
%             end
            this.Ma = bestMa;

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


