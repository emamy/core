classdef KernelInterpol < KerMorObject & IKernelCoeffComp
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
    % @change{0,7,dw,2014-01-24} Removed the preconditioning stuff and LU
    % factorization. Instead included the Newton basis computation method
    % from \cite PS11
    %
    % @change{0,7,dw,2013-01-23} Re-added the LU decomposition stuff to this class from
    % MemoryMatrix, as the class has been removed.
    %
    % @new{0,6,dw,2012-01-23} Included preconditioning techniques for kernel interpolation
    % according to \cite S08.
    %
    % @change{0,5,dw,2011-09-12} Set the UseLU flag to true per default.
    % Using FileMatrix instances now, along with flags of whether to
    % successively build the inverse, too.
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
    
    properties(SetObservable)
        % Flag that indicates whether to apply the Newton basis for stable
        % interpolation \cite PS11 .
        %
        % @propclass{optional} This can greatly improve numerical stability
        % of computation, but possibly takes more time to compute.
        %
        % @default true
        UseNewtonBasis = true;
    end
    
    properties(Access=private)
        % The kernel matrix K to use within interpolation.
        %
        % @propclass{data} Required for any interpolation computation.
        %
        % @type data.FileMatrix
        K;
    end
    
    properties(Access=private)
        % The cholesky decomposition matrix for Newton basis case
        kexp;
    end
    
    methods
        
        function this = KernelInterpol
            this = this@KerMorObject;
            this.MultiTargetComputation = true;
            this.registerProps('UseNewtonBasis');
        end
        
        function copy = clone(this)
            copy = general.interpolation.KernelInterpol;
            copy = clone@IKernelCoeffComp(this, copy);
            copy.UseNewtonBasis = this.UseNewtonBasis;
            copy.K = this.K;
            copy.kexp = this.kexp;
        end
        
        function [a, base, used] = interpolate(this, fxi)
            % Computes the kernel expansion coefficients `\alpha_i`.
            %
            % Parameters:
            % fxi: The real function value samples at centers `x_i`
            %
            % Return values:
            % a: The coefficient vector `\alpha` @type matrix<double>
            % base: The newton basis values to set for
            % kernels.KernelExpansion.Base @type matrix<double>
            N = size(fxi,2);
            if this.UseNewtonBasis   
                NV = zeros(N,N);
                c = zeros(size(fxi,1),N);
                fxin = Norm.L2(fxi);
                sumNsq = zeros(1,N);
                fresidual = fxi;
                used = zeros(1,N);
                for m = 1:N
                    res = sum(fresidual.^2,1);
                    [~, maxidx] = max(res);
                    col = this.kexp.getKernelMatrixColumn(maxidx);
                    tN = col - NV(:,1:m-1)*NV(maxidx,1:m-1)';
                    tNnorm = sqrt(col(maxidx)-sumNsq(maxidx));
                    if max(res./fxin) < 1e-10 || ~isreal(tNnorm) || tNnorm <= 0
                        used = used(1:m-1);
                        break;
                    end
                    NV(:,m) = tN/tNnorm;
                    c(:,m) = fresidual(:,maxidx)./tNnorm;                    
                    fresidual = fresidual - c(:,m)*(NV(:,m))';
                    sumNsq = sumNsq + (NV(:,m).^2)';
                    used(m) = maxidx;
                end
                base = eye(N);
                base(:,1:length(used)) = NV(:,used);
                %a = zeros(size(c));
                %a(:,1:length(used)) = c(:,used);
                a = c(:,used);
            else
                base = 1;
                a = (this.K\fxi')';
                used = 1:N;
            end
        end
        
        %% IKernelCoeffComp interface members
        function init(this, kexp)
            % Initializes the interpolation
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            this.K = [];
            this.kexp = [];
            if this.UseNewtonBasis
                this.kexp = kexp;
            else
                this.K = kexp.getKernelMatrix;
            end
        end
        
        function [ci, svidx, sf] = computeKernelCoefficients(this, fxi, ~)
            % Implementation of the kernels.ICoeffComp interface
            %
            % Parameters:
            % fxi: The target values `\vf(\vx_i)` as row vector @type rowvec<double>
            %
            % Return values:
            % ci: The coefficients `c_i` as row vector @type rowvec
            % svidx: The support vector indices of all elements of `c_i`
            % that regarded to be support vectors. @type integer
            
            [ci, base, svidx] = this.interpolate(fxi);
            % Here: Transform the coefficients back to direct translate
            % basis as the Base cannot be set here and might be different
            % for each component function
            if this.UseNewtonBasis
                ci = ci / base;
            end
            sf = StopFlag.SUCCESS;
        end
    end
    
    methods(Static)
        
        function test_KernelInterpolation
            % Performs a test of this class
            
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            
%             x = -5:.05:5;
%             fx = sinc(x);
            
            x = 0:.05:10;
            fx = x.^2.*sin(2*x);
            
            n = length(x);
            samp = 1:15:n;
            
            kexp = kernels.KernelExpansion;
            
            kexp.Kernel = kernels.GaussKernel(.5);
            internal_test(x,fx,false);
            internal_test(x,fx,true);
            
            internal_test(x,ones(size(x))*5,false);
            internal_test(x,ones(size(x))*5,true);
            
            kexp.Kernel = kernels.InvMultiquadrics(-1.4,2);
            internal_test(x,fx,true);
            internal_test(x,fx,false);
            
            kexp.Kernel = kernels.InvMultiquadrics(-4,5);
            internal_test(x,fx,true);
            internal_test(x,fx,false);
            
            pm.done;
            
            function internal_test(x,fx,lu)
                xi = x(samp);
                fxi = fx(samp);
                
                ki = general.interpolation.KernelInterpol;
                
                kexp.Centers.xi = xi;
                ki.UseNewtonBasis = lu;
                ki.init(kexp);
                
                h = pm.nextPlot('xx',sprintf('Interpolation test with kernel %s, newton=%d',class(kexp.Kernel),lu));
                plot(h,x,fx,'r');
                
                [kexp.Ma, kexp.Base] = ki.interpolate(fxi);
                
                fsvr = kexp.evaluate(x);
                
                hold(h,'on');
                % Plot approximated function
                plot(h,x,fsvr,'b--');
                %skipped = setdiff(1:length(x),svidx);
                plot(h,xi,fxi,'.r');
            end
        end
        
        function test_KernelInterpolationError
            realprec = 20;
            range = [-4 4];
            
            pm = PlotManager(false,3,3);
            pm.LeaveOpen = true;
            
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernels.GaussKernel(1);
            testfun = @(x)sinc(x) * 2 .*sin(3*x);
            %             kernel = kernels.InvMultiquadrics(-1,1);
            %             testfun = @(x)x.^3+2*x.^2-3*x;
            
            ki = general.interpolation.KernelInterpol;
            
            hsteps = [.5:.25:2 2:2:10];
            err = zeros(size(hsteps));
            for n = 1:length(hsteps)
                h = 1/hsteps(n);
                xi = range(1):h:range(2);
                xi = xi + rand(size(xi))*max(xi)/100;
                fxi = testfun(xi);
                
                x = range(1):h/realprec:range(2);
                fx = testfun(x);
                
                kexp.Centers.xi = xi;
                ki.init(kexp);
                [kexp.Ma, kexp.Base] = ki.interpolate(fxi);
                fint = kexp.evaluate(x);
                
                ax = pm.nextPlot(['step' n],sprintf('h=%g',h),'x','fx');
                plot(ax,xi,fxi,'r',x,fx,'r--',x,fint,'b');
                err(n) = max(abs(fx-fint));
            end
            ax = pm.nextPlot('errors','Error development over h sizes','h','errors');
            h = 1./hsteps;
            plot(ax,h,exp(log(h)./h),'b',h,exp(-1./h),'b--',h,err,'r');
            legend(ax,'exp(log(h)/h)','exp(-1/h)','||f-sfx||_\infty');
        end
    end
end

