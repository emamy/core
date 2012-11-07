classdef KernelInterpol < KerMorObject & approx.algorithms.IKernelCoeffComp
    % Provides kernel interpolation.
    %
    % The basic interpolation form is 
    % `` f(x) = \sum\limits_{i=1}^N \alpha_i \Phi(x,x_i)``
    % Interpolation finds coefficients such that 
    % `fx_i = RHS(x_i)` for `i=1\ldots N`.
    %
    % There is also a zero-function threshold `10*eps`. If all
    % `fx_i-\beta` values are below that a constant function is
    % assumed.
    %
    % The preconditioning technique is implemented after \cite S08.
    %
    % @author Daniel Wirtz @date 01.04.2010
    %
    % @new{0,6,dw,2012-01-23} Included preconditioning techniques for kernel interpolation
    % according to \cite S08.
    %
    % @change{0,5,dw,2011-09-12} Set the UseLU flag to true per default.
    % Using IKernelMatrix instances now, along with flags of whether to
    % successively build the inverse, too. Moved the UseLU property to
    % MemoryKernelMatrix.
    %
    % @change{0,4,dw,2011-05-03} Removed the artificial offset term `b` from the interpolation
    % process (no longer used in kernel expansions)
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @change{0,2,dw,2011-03-21} 
    % - Added the general.interpolation.KernelInterpol.UseLU property. With
    % this subsequent calls to interpolate using the same kernel matrix is
    % more efficient.
    % - Updated the documentation
    
    properties(Dependent)
        % The kernel matrix K to use within interpolation.
        %
        % @propclass{data} Required for any interpolation computation.
        %
        % @type data.IKernelMatrix
        K;
    end
    
    properties(Access=private)
        fK;
        
        % Preconditioning matrix
        P;
    end
    
    methods
        
        function this = KernelInterpol
            this = this@KerMorObject;
            this.MultiTargetComputation = true;
            this.registerProps('K');
        end
        
        function a = interpolate(this, fxi)
            % Computes the kernel expansion coefficients `\alpha_i`.
            %
            % Parameters:
            % fxi: The real function value samples at centers `x_i`
            %
            % Return values:
            % a: The coefficient vector `\alpha`

            if isa(this.fK, 'data.MemoryKernelMatrix')
                if this.fK.UseLU
                    if ~isempty(this.P)
                        warning('KerMor:interpol','Preconditioning matrix available, not applying as UseLU is set!');
                    end
                    a = this.fK.U\(this.fK.L\fxi');
                    return;
                elseif this.fK.BuildInverse
%                     P = this.fK.Kinv * this.fK.K;
%                     pfxi = this.fK.Kinv*fxi';
%                     a = P\pfxi;
                    if ~isempty(this.P)
                        warning('KerMor:interpol','Preconditioning matrix available, not applying as BuildInverse is set!');
                    end
                    a = this.fK.Kinv * fxi';
                    return;
                end
            end
            %warning('MATLAB:nearlySingularMatrix','off');
            if ~isempty(this.P)
                % fK is already preconditioned, see init(..)
                a = this.fK\(this.P*fxi');
            else
                a = this.fK\fxi';
            end
        end
        
        function set.K(this, value)
            if ismatrix(value) && isa(value,'double')
                value = data.MemoryKernelMatrix(value);
            end
            if ~isa(value, 'data.IKernelMatrix')
                error('K must be a data.IKernelMatrix or a double matrix.');
            end
            this.fK = value;
        end
        
        function K = get.K(this)
            K = this.fK;
        end
        
        %% approx.algorithms.IKernelCoeffComp interface members
        function init(this, K, kexp)
            % Initializes the interpolation
            %
            % Parameters:
            % K: The kernel matrix @type data.MemoryKernelMatrix
            % kexp: If given, some preconditioning can be performed usign the kernel expansions
            % information. @type kernels.KernelExpansion @default []
            this.K = K;
            this.P = [];
            if nargin > 2
                if ~isa(this.fK,'data.MemoryKernelMatrix')
                    error('Preconditioning with other than MemoryKernelMatrices not yet implemented.');
                end
                if isa(kexp.Kernel,'kernels.ARBFKernel') && kexp.Kernel.epsilon < .01
                    [this.P, k2] = this.getPreconditioner(kexp.Kernel, kexp.Centers.xi);
                    % Overwrite current matrix with preconditioned one
                    oldK = this.fK.K;
                    this.fK = data.MemoryKernelMatrix(this.P*this.fK.K);
                    x = kexp.Centers.xi;
                    [m M] = general.Utils.getBoundingBox(x);
                    c1 = cond(oldK);
                    c2 = cond(this.fK.K);
                    pl = '-';
                    if c2<c1
                        pl = '+';
                    end
                    if KerMor.App.Verbose > 3
                        fprintf('Cond(K)=%e, Cond(P*K)=%e, size(x)=[%d %d], %s! eps=%e, xdiag=%e, k2=%d\n',...
                            c1,c2,size(x,1),size(x,2),pl,kexp.Kernel.epsilon,norm(M-m),k2);
                    end
                end
            end
        end
        
        function [ci, svidx] = computeKernelCoefficients(this, fxi, ~)
            % Implementation of the kernels.ICoeffComp interface
            %
            % Parameters:
            % fxi: The target values `\vf(\vx_i)` as row vector @type rowvec<double>
            %
            % Return values:
            % ci: The coefficients `c_i` as row vector @type rowvec
            % svidx: The support vector indices of all elements of `c_i`
            % that regarded to be support vectors. @type integer
            
            % Transform to row vector
            ci = this.interpolate(fxi)';
            svidx = 1:size(ci,2);
        end
    end
    
    
    
    methods(Static)
        function [P, k2] = getPreconditioner(k, x)
            % Computes the preconditioning matrix `D_\ep M` from the work \cite S08 given a kernel
            % expansion and the centers at which the function values are to be interpolated.
            %
            % Parameters:
            % k: A radial kernel `\Phi`. @type kernels.ARBFKernel
            % x: The centers `x_1 \ldots x_N \in\R^n` at which to interpolate using the kernel `\Phi`
            %
            % Return values:
            % P: The preconditioning matrix `P\in\R^{n\times n}`
            % k2: The number `k_2` of \cite S08
            
            if ~isa(k,'kernels.ARBFKernel')
                error('The kernel used must be a kernels.ARBFKernel subclass');
            end
            
            N = size(x,2);
            M = [];
            I_N = eye(N);
            pivi = 1;
            cols = 1;
            L = I_N; P = I_N;
            tk = zeros(1,N);
            mi = general.MonomialIterator(size(x,1));
            
            while pivi <= N
                % Compute next monomial
                
                % Stragegy 1: Random
                %                     deg = ceil(randn(1)*N);
                %                     alpha = mi.getRandMonomial(deg);
                
                % Stragegy 2: Ordered list
                if pivi == 1
                    alpha = mi.getNullMonomial;
                    deg = 0;
                else
                    [alpha, deg] = mi.nextMonomial;
                end
                
                % Compute new moment matrix column
                newMcol = prod(x .^ repmat(alpha,1,N),1)';
                M = [M newMcol];%#ok
                
                % apply previous changes to new column
                newMcol = L*P*newMcol;
                
                % Select pivot element candidate indices in current column
                sel = pivi:N;
                
                % Strategy one: Use maximum pivoting
                %                     [v, maxidx] = max(abs(newMcol(sel)));
                %                     permidx = sel(maxidx);
                
                % Strategy two: Only find first nonzero-row and use it
                permidx = sel(find(abs(newMcol(sel)) > sqrt(eps),1));
                if ~isempty(permidx)
                    v = abs(newMcol(permidx));
                else
                    v = 0;
                end
                
                % step one column ahead if current column is already annihilated
                if v < sqrt(eps)
                    %v
                    cols = cols+1;
                    continue;
                end
                
                % get new permutation matrix according to pivot
                Pn = getPermMat(pivi, permidx, N);
                % swap columns
                newMcol = Pn*newMcol;
                % compose L_k
                Ln = I_N;
                l_k = -newMcol(pivi+1:N)/newMcol(pivi);
                Ln(pivi+1:N,pivi) = l_k;
                
                % keep record of the column indices at which the next linear independent
                % monomial was added
                tk(pivi) = deg;
                
                L = Ln*Pn*L*Pn;
                P = Pn*P;
                
                pivi = pivi+1;
                cols = cols+1;
            end
            
            %U = L*P*M;
            D = diag(k.epsilon.^-tk);
            %PM = P; % return accum. permutation matrix
            P = D * L* P; % compute preconditioning matrixl
            k2 = tk(end);
            
            function P = getPermMat(i,j,n)
                P = eye(n);
                P(i,i) = 0; P(j,j) = 0;
                P(j,i) = 1; P(i,j) = 1;
            end
        end   
        
% ------------ USED IN PAPER WH10 synth. model tests
%         function test_KernelInterpolation2()
%             dim = 3;
%             x = repmat(0:pi/10:50,dim,1);
%             fx = -sin(.5*x).*x*.2;
%             fx = fx(1,:);
%             
%             n = size(x,2);
%             samp = 11:20:n;
%             xi = x(:,samp);
%             fxi = fx(samp);
%             
%             kernel = kernels.GaussKernel(1);
%             ki = general.interpolation.KernelInterpol;
%             figure(1);
%             plot(x(1,:),fx,'r');
%             
%             g = 1:.1:3*pi^2;
%             minidx = 0; mina = Inf; mintest = Inf;
%             for idx = 1:length(g)
%                 kernel.Gamma = g(idx);
%                 ki.K = kernel.evaluate(xi,xi);
%                 
%                 [a,b] = ki.interpolate(fxi);
%                 
%                 svfun = @(x)a'*kernel.evaluate(xi,x)+b;
% 
%                 fsvr = svfun(x);
%                 
%                 if norm(fx-fsvr) < mina
%                     mina = norm(fx-fsvr);
%                     minidx=idx;
%                     minfsvr = fsvr;
%                 end
%                 if norm(fx-fsvr)*norm(a) < mintest
%                     minidxtest=idx;
%                     minfsvrtest = fsvr;
%                     minav = a;
%                     minb = b;
%                     mintest = norm(fx-fsvr)*norm(a);
%                 end
% 
% %                 hold on;
% % 
% %                 % Plot approximated function
% %                 plot(x,fsvr,'b--');
% %                 %skipped = setdiff(1:length(x),svidx);
% %                 plot(xi,fxi,'.r');
% % 
% %                 hold off;
%                 
%                 fprintf('g=%f, anorm=%f, diff=%f\n',g(idx),norm(a),norm(fx-fsvr));
%             end
%             hold on;
%             plot(x,minfsvr,'b--',x,minfsvrtest,'g--');
%             %skipped = setdiff(1:length(x),svidx);
%             plot(xi,fxi,'.r');
%             hold off;
%             disp(g(minidx));
%             disp(g(minidxtest));
%             minav
%             minb
%         end
        
        function test_KernelInterpolation()
            % Performs a test of this class
            
%             x = -5:.05:5;
%             fx = sinc(x);
            
            x = 0:.05:10;
            fx = x.^2.*sin(2*x);
            
            n = length(x);
            samp = 1:15:n;
            
            k = kernels.GaussKernel(.5);
            internal_test(x,fx,samp,k,10);
            
            internal_test(x,ones(size(x))*5,samp,k,13);
            
            k = kernels.InvMultiquadrics(-1.4,2);
            internal_test(x,fx,samp,k,11);
            
            k = kernels.InvMultiquadrics(-4,5);
            internal_test(x,fx,samp,k,12);
            
            function internal_test(x,fx,samp,kernel,fignr)
                xi = x(samp);
                fxi = fx(samp);
                
                ki = general.interpolation.KernelInterpol;
                
                ki.K = kernel.evaluate(xi,xi);
                
                figure(fignr);
                plot(x,fx,'r');
                
                a = ki.interpolate(fxi);
                
                svfun = @(x)a'*kernel.evaluate(xi,x);
                
                fsvr = svfun(x);
                
                hold on;
                
                % Plot approximated function
                plot(x,fsvr,'b--');
                %skipped = setdiff(1:length(x),svidx);
                plot(xi,fxi,'.r');
                
                hold off;
            end
        end
        
        function test_KernelInterpolationError
            realprec = 20;
            range = [-4 4];
            
            kernel = kernels.GaussKernel(2);
            testfun = @(x)sinc(x) * 2 .*sin(3*x);
            %             kernel = kernels.InvMultiquadrics(-1,1);
            %             testfun = @(x)x.^3+2*x.^2-3*x;
            
            ki = general.interpolation.KernelInterpol;
            
            hsteps = [.5:.25:2 2:2:10];
            err = zeros(size(hsteps));
            for n = 1:length(hsteps)
                h = 1/hsteps(n);
                xi = range(1):h:range(2);
                fxi = testfun(xi);
                
                x = range(1):h/realprec:range(2);
                fx = testfun(x);
                
                ki.K = kernel.evaluate(xi,xi);
                a = ki.interpolate(fxi);
                fint = a'*kernel.evaluate(xi,x);
                
                err(n) = max(abs(fx-fint));
                
                %                 figure(n+1);
                %                 plot(xi,fxi,'r.',x,fx,'r',x,fint,'b--',x,abs(fx-fint),'g');
                %                 title(sprintf('h=%f, ||f-sfx||=%f',h,err(n)));
            end
            figure(1);
            h = 1./hsteps;
            plot(h,exp(log(h)./h),'b',h,exp(-1./h),'b--',h,err,'r');
            legend('exp(log(h)/h)','exp(-1/h)','||f-sfx||_\infty');
            xlabel('h');
        end
    end
end

