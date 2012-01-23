classdef GaussKernel < kernels.BellFunction
    % Radial Basis Function Kernel
    %   
    % Uses the notation
    % ``\Phi(x,y) = e^{\frac{||x-y||^2}{\gamma^2}},``
    % so be careful with the `\gamma` constant.
    %
    % @author Daniel Wirtz @date 11.03.2011
    %
    % @change{0,5,dw,2011-10-16} Exported the evaluate function to
    % kernels.ARBFKernel, but re-implemented the customized
    % evaluate function as the norm squared is already computed fast and
    % first taking the square root and then squaring again would introduce
    % unecessary overhead.
    %
    % @change{0,3,dw,2011-04-26} Fixed the x0 computation for the new Gamma property version; so far
    % a square was missing, rendering the LipschitzEstimator modified newton method useless.
    %
    % @change{0,3,dw,2011-04-15} Changed the Gamma property to compute into the kernel evaluation
    % being squared instead of linear. This way the Gamma becomes a more geometrical meaning
    %
    % @change{0,2,dw,2011-03-11} Added new speed tests for one and two
    % argument calls to 'evaluate'. The tests are run 'iter' times and the
    % mean value is plotted to the output. (See @ref
    % kernels.GaussKernel.test_GaussMexSpeedTest1Arg and @ref
    % kernels.GaussKernel.test_GaussMexSpeedTest2Arg)
    
    properties(SetObservable)
        % Univariate scaling
        %
        % @propclass{critical} Greatly influences the kernels behaviour.
        %
        % @type double @default 1
        Gamma = 1;
    end    
    
    methods
        function this = GaussKernel(Gamma)
            % Creates a new GaussKernel
            %
            % Parameters:
            % Gamma: The Gamma property to use. @default 1 @type double
            this = this@kernels.BellFunction;
            this.registerProps('Gamma');%,'Sigma'
            
            if nargin == 1
                this.Gamma = Gamma;
            else
                % Set r0 according to whatever default value for Gamma is currently set
                this.r0 = sqrt(this.Gamma^2/2);
            end
        end
        
        function K = evaluate(this, x, y)
            % Evaluates the gaussian.
            %
            % If `y_j` is set, the dimensions of `x_i` and `y_j` must be equal for all `i,j`.
            %
            % Parameters:
            % x: First set `x_i \in \R^d` of `n` vectors @type matrix<double>
            % y: Second set `y_j \in \R^d` of `m` vectors. If y is empty `y_i = x_i` and `n=m`
            % is assumed. @type matrix<double>
            %
            % Return values:
            % K: An evaluation matrix `K \in \R^{n\times m}` of the evaluated gaussians with
            % entries `K_{i,j} = e^{-\norm{x_i-y_j}{G}^2/\gamma^2}`.
            K = exp(-this.getSqDiffNorm(x, y)/this.Gamma^2);
        end
                
        function Nablax = getNabla(this, x, y)
            % Method for first derivative evaluation
            if size(x,2) > 1 && size(y,2) > 1
                error('One argument must be a vector.');
            end
            if ~isempty(this.P)
                error('Not yet implemented correctly.');
%                 xl = x(this.P,:);
%                 yl = y(this.P,:);
%             else
%                 xl = x; yl = y;
            end
            hlp = bsxfun(@minus,x,y);
            hlp = -2*hlp/this.Gamma^2;
            Nablax = bsxfun(@times,hlp,this.evaluate(x, y));
        end 
        
%         function K = evaluateIntel(this, x, varargin)
%             % Experimental function that automatically calls the mex openmp
%             % implementation code if the vectors are small enough.
%             %
%             % @todo write c code more efficient (use blas/lapack?)
%             
%             % Evaluate MEX function if sizes are small!
%             if numel(x) < 500000
%                  K = this.evaluateMex(x,varargin{:});
%                  return;
%             end
%             
%             n1sq = sum(x.^2,1);
%             n1 = size(x,2);
% 
%             if nargin == 2;
%                 n2sq = n1sq;
%                 n2 = n1;
%                 y = x;
%             else
%                 y = varargin{1};
%                 n2sq = sum(y.^2,1);
%                 n2 = size(y,2);
%             end;
%             K = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq - 2*x'*y;
%             K(K<0) = 0;
%             K = exp(-K/this.Gamma^2);
%         end
                
        function dx = evaluateD1(this, r)
            % Method for first derivative evaluation
            dx = -2*r/this.Gamma^2 .* exp(-r.^2/this.Gamma^2);
        end
        
        function ddx = evaluateD2(this, r)
            % Method for second derivative evaluation
            ddx = (2/this.Gamma^2) * (2*r.^2/this.Gamma^2-1) .* exp(-r.^2/this.Gamma^2);
        end
        
        function phi = evaluateScalar(this, r)
            % Implements the required method from the IRotationInvariant
            % interface
            phi = exp(-r.^2/this.Gamma^2);
        end
        
        function set.Gamma(this, value)
            % @todo check why penalty factor was set here to 1/value! ?!?
            if ~isreal(value) || ~isscalar(value) || value <= 0
                error('Only positive scalar values allowed for Gamma.');
            end
            this.Gamma = value;
            
            % Update the dilation parameter of ARBFKernel
            this.epsilon = 1/value;
            
            % Adjust the BellFunctions' x0 value
            this.r0 = value/sqrt(2);
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
            g = dist/sqrt(-log(ep));
            this.Gamma = g;
            
            if KerMor.App.Verbose > 3
                fprintf('Setting Gamma = %12.20f\n',g);
            end
        end
    end
    
    methods(Static)
        
        function res = test_InterpolGamma            
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

