classdef GaussKernel < kernels.BellFunction
    % Radial Basis Function Kernel
    %   
    % Uses the notation
    % ``\Phi(x,y) = e^{\frac{||x-y||^2}{\gamma}},``
    % so be careful with the `\gamma` constant.
    %
    % @author Daniel Wirtz @date 11.03.2011
    %
    % @change{0,2,dw,2011-03-11} Added new speed tests for one and two
    % argument calls to 'evaluate'. The tests are run 'iter' times and the
    % mean value is plotted to the output. (See @ref
    % kernels.GaussKernel.test_GaussMexSpeedTest1Arg and @ref
    % kernels.GaussKernel.test_GaussMexSpeedTest2Arg)
    
    properties
        Gamma = 1;
    end
    
    methods
        function this = GaussKernel(Gamma)
            if nargin == 1
                this.Gamma = Gamma;
            end
        end
        
        function K = evaluate(this, x, y)
            n1sq = sum(x.^2,1);
            n1 = size(x,2);

            if nargin == 2;
                n2sq = n1sq;
                n2 = n1;
                y = x;
            else
                n2sq = sum(y.^2,1);
                n2 = size(y,2);
            end;
            K = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq - 2*x'*y;
            K(K<0) = 0;
            K = exp(-K/this.Gamma);
        end
        
        function K = evaluateIntel(this, x, varargin)
            % Experimental function that automatically calls the mex openmp
            % implementation code if the vectors are small enough.
            %
            % @todo write c code more efficient (use blas/lapack?)
            
            % Evaluate MEX function if sizes are small!
            if numel(x) < 500000
                 K = this.evaluateMex(x,varargin{:});
                 return;
            end
            
            n1sq = sum(x.^2,1);
            n1 = size(x,2);

            if nargin == 2;
                n2sq = n1sq;
                n2 = n1;
                y = x;
            else
                y = varargin{1};
                n2sq = sum(y.^2,1);
                n2 = size(y,2);
            end;
            K = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq - 2*x'*y;
            K(K<0) = 0;
            K = exp(-K/this.Gamma);
        end
                
        function dx = evaluateD1(this, x)
            % Method for first derivative evaluation
            dx = -2*x/this.Gamma .* exp(-x.^2/this.Gamma);
        end
        
        function ddx = evaluateD2(this, x)
            % Method for second derivative evaluation
            ddx = (2/this.Gamma) * (2*x.^2/this.Gamma-1) .* exp(-x.^2/this.Gamma);
        end
        
        function Kx = evaluateScalar(this, x)
            % Implements the required method from the IRotationInvariant
            % interface
            Kx = exp(-x.^2/this.Gamma);
        end
        
        function set.Gamma(this, value)
            if ~isposrealscalar(value)
                error('Only positive scalar values allowed for Gamma.');
            end
            this.Gamma = value;
            % Adjust the BellFunctions' x0 value
            this.x0 = sqrt(value/2); %#ok
            this.PenaltyFactor = 1/value; %#ok
        end
        
        function g = setGammaForDistance(this, dist, ep)
            % Computes the `\gamma` value for which the Gaussian is smaller
            % than `\epsilon` in a distance of dist, i.e.
            % ``e^{-\frac{d^2}{\gamma}) < \epsilon``
            % Returns the computed value AND sets the kernel's Gamma
            % property to this value.
            %
            % Parameters:
            % dist: The target distance at which the gaussian is smaller
            % than ep
            % ep: The `\epsilon` value. If not given, `\epsilon`=eps
            % (machine precision) is assumed.
            %
            % Return values:
            % g: The computed gamma
            if nargin == 2
                ep = eps;
            end
            g = -(dist^2)/log(ep);
            this.Gamma = g;
            
            if KerMor.App.Verbose > 0
                fprintf('Setting Gamma = %f\n',g);
            end
        end
    end
    
    methods(Static)
        
        function res = test_InterpolGamma
            
% Algorithmus
% evaluate f over whole validationdata (subset projdata) and store
% select initial training set from projdata
%   find optimal kernel config (for all dims!)
%   compute position of max error (normdiff f-evals)
%   add pos to training set
% while numcenters < maxnum || err < tol
            
            ki = general.interpolation.KernelInterpol;
            %ki.UseLU = true;
            k = kernels.GaussKernel;
            dx = .2;
            x = -3:dx:3;
            %xfine = -3:dx/2:3;
            fx = sin(x*pi); %+rand(size(x))
            plot(x,fx);
            epsteps = 0.05:.05:.95;
            dlog = zeros(3,length(epsteps));
            for epidx=1:length(epsteps)
                ep = epsteps(epidx);
                k.setGammaForDistance(dx,ep);
                for idx = 1:length(x)
                    x2 = x;
                    x2(idx) = [];
                    fx2 = fx;
                    fx2(idx) = [];
                    
                    ki.K = k.evaluate(x2);
                    [a,b] = ki.interpolate(fx2);
                    
                    fxi = a'*k.evaluate(x2,x) + b;
                    diff(idx) = abs(fx(idx) - fxi(idx));
                end
                dlog(1,epidx) = ep;
                dlog(2,epidx) = min(diff);
                dlog(3,epidx) = max(diff);
            end            
            disp(dlog);
            [val, idx] = min(dlog(2,:));
            fprintf('Min distance: %f at ep=%f\n',dlog(2,idx),dlog(1,idx));
        end
        
        function res = test_GaussMexSpeedTest1Arg(sx,sy,iter)
            if nargin < 3
                iter = 50;
                if nargin < 2
                    sy = 100;
                    if nargin < 1
                        sx = 5000;
                    end
                end
            end
            k = kernels.GaussKernel(1);            
            x = rand(sx,sy);
            
            fprintf('One argument speed test with sx=%d, sy=%d and %d iterations\n',sx,sy,iter);
            tmex = zeros(1,iter); tmex2 = zeros(1,iter);
            tmexp = zeros(1,iter); tmat = zeros(1,iter);
            fprintf('Iteration ');
            for i=1:iter
                fprintf('%d ',i);
                
                t = tic;
                Kmex = k.dontuse_evaluate(x);
                %Kmex = k.evaluateIntel(x);
                tmex(i) = toc(t);

                t = tic;
                Kmex2 = k.dontuse_evaluateDirect(x);
                tmex2(i) = toc(t);

                t = tic;
                KmexP = k.evaluateMex(x);
                tmexp(i) = toc(t);

                t = tic;
                K = k.evaluate(x);
                tmat(i) = toc(t);
            
            end
            fprintf('done!\n');
            fprintf(['1: %1.5fs - Mex straight\n2: %1.5fs - Mex time opt\n'...
                '3: %1.5fs - Mex time opt openmp\n4: %1.5fs - Matlab time\n'...
                'Difference norms: 1-4=%1.5f, 2-4=%1.5f, 3-4=%1.5f\n'],...
                mean(tmex2),mean(tmex),mean(tmexp),mean(tmat),...
                norm(Kmex2-K),norm(Kmex-K),norm(KmexP-K));
            
            res = true;
        end
        
        function res = test_GaussMexSpeedTest2Arg(sx,sy1,sy2,iter)
            if nargin < 4
                iter = 40;
                if nargin < 3
                    sy2 = 100;
                    if nargin < 2
                        sy1 = 100;
                        if nargin < 1
                            sx = 500;
                        end
                    end
                end
            end
            
            k = kernels.GaussKernel(1);            
            x = rand(sx,sy1);
            y = rand(sx,sy2);
            
            fprintf('Two argument speed test with sx=%d, sy1=%d, sy2=%d and %d iterations\n',sx,sy1,sy2,iter);
            tmex = zeros(1,iter); tmex2 = zeros(1,iter);
            tmexp = zeros(1,iter); tmat = zeros(1,iter);
            fprintf('Iteration ');
            for i=1:iter
                fprintf('%d ',i);
                
                t = tic;
                Kmex = k.dontuse_evaluate(x,y);
                %Kmex = k.evaluateIntel(x,y);
                tmex(i) = toc(t);

                t = tic;
                Kmex2 = k.dontuse_evaluateDirect(x,y);
                tmex2(i) = toc(t);

                t = tic;
                KmexP = k.evaluateMex(x,y);
                tmexp(i) = toc(t);

                t = tic;
                K = k.evaluate(x,y);
                tmat(i) = toc(t);
            
            end
            fprintf('done!\n');
            fprintf(['1: %1.5fs - Mex straight\n2: %1.5fs - Mex time opt\n'...
                '3: %1.5fs - Mex time opt openmp\n4: %1.5fs - Matlab time\n'...
                'Difference norms: 1-4=%1.5f, 2-4=%1.5f, 3-4=%1.5f\n'],...
                mean(tmex2),mean(tmex),mean(tmexp),mean(tmat),...
                norm(Kmex2-K),norm(Kmex-K),norm(KmexP-K));
            
            res = true;
        end
    end
    
end

