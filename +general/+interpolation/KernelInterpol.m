classdef KernelInterpol < KerMorObject & IKernelCoeffComp
    % Provides kernel interpolation.
    %
    % The basic interpolation form is 
    % `` f(x) = \sum\limits_{i=1}^N \alpha_i \K(x,x_i)``
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
    % from \cite PS11 and controllable interpolation precision via RelTol
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
        % interpolation \cite PS11 . Using this interpolation scheme also
        % allows to control the interpolation precision via setting RelTol.
        %
        % @propclass{optional} This can greatly improve numerical stability
        % of computation, but possibly takes more time to compute.
        %
        % @type logical @default true
        %
        % See also: RelTol
        UseNewtonBasis = true;
        
        % Maximum relative error tolerance (L2 in function dimension) over
        % the given training data. Only relevant when UseNewtonBasis =
        % true. This setting causes the iterative Newton algorithm to
        % stop as soon as the maximum relative error over the training data
        % is less than RelTol, e.g. ensures an interpolation on the
        % training data up to the specified precision.
        %
        % This is useful in particular, if the training data is
        % highly/almost redundant and the direct inversion is too
        % numericall ill conditioned.
        %
        % @type double @default 1e-13
        %
        % @propclass{important} Specifying a higher tolerance results in
        % sparser but less precise approximations
        %
        % See also: UseNewtonBasis
        RelTol = 1e-13;
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
        
        function [ai, nbase, usedcenteridx] = interpolate(this, fxi)
            % Computes the kernel expansion coefficients `\alpha_i`.
            %
            % Parameters:
            % fxi: The real function value samples at centers `x_i`
            %
            % Return values:
            % ai: The coefficient vector `\alpha` @type matrix<double>
            % nbase: The newton basis values to set for
            % kernels.KernelExpansion.Base @type matrix<double>
            % usedcenteridx: The indices of the used centers.
            N = size(fxi,2);
            if this.UseNewtonBasis   
                NV = zeros(N,N);
                c = zeros(size(fxi,1),N);
                fxin = Norm.L2(fxi);
                sumNsq = zeros(1,N);
                fresidual = fxi;
                usedcenteridx = zeros(1,N);
                for m = 1:N
                    res = sum(fresidual.^2,1);
                    [~, maxidx] = max(res);
                    col = this.kexp.getKernelMatrixColumn(maxidx);
                    tN = col - NV(:,1:m-1)*NV(maxidx,1:m-1)';
                    tNnorm = sqrt(col(maxidx)-sumNsq(maxidx));
                    if max(res./fxin) < this.RelTol || ~isreal(tNnorm) || tNnorm <= 0
                        m = m-1;%#ok
                        break;
                    end
                    NV(:,m) = tN/tNnorm;
                    c(:,m) = fresidual(:,maxidx)./tNnorm;                    
                    fresidual = fresidual - c(:,m)*(NV(:,m))';
                    sumNsq = sumNsq + (NV(:,m).^2)';
                    usedcenteridx(m) = maxidx;
                end
                usedcenteridx = usedcenteridx(1:m);
                nbase = NV(:,1:m);
                ai = c(:,1:m);
                
                usedcenteridx = sort(usedcenteridx,'ascend');
                nbase = nbase(usedcenteridx,:);
                if length(usedcenteridx) < N && nargout < 3
                    error('Less basis translates required than given, but their indices are not considered as output. Please receive the "used" output and select the indicated subset as centers.');
                end
            else
                nbase = 1;
                ai = (this.K\fxi')';
                usedcenteridx = 1:N;
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
            
            ke = kernels.KernelExpansion;
            
            ke.Kernel = kernels.GaussKernel(.5);
            internal_test(x,fx,false);
            internal_test(x,fx,true);
            
            internal_test(x,ones(size(x))*5,false);
            internal_test(x,ones(size(x))*5,true);
            
            ke.Kernel = kernels.InvMultiquadrics(-1.4,2);
            internal_test(x,fx,true);
            internal_test(x,fx,false);
            
            ke.Kernel = kernels.InvMultiquadrics(-4,5);
            internal_test(x,fx,true);
            internal_test(x,fx,false);
            
            % 2D test
            fx = [fx; x.^3.*sin(.4*x)];
            internal_test(x,fx,true);
            internal_test(x,fx,false);
            
            pm.done;
            
            function internal_test(x,fx,lu)
                xi = x(samp);
                fxi = fx(:,samp);
                
                ki = general.interpolation.KernelInterpol;
                
                ke.Centers.xi = xi;
                ki.UseNewtonBasis = lu;
                ki.init(ke);
                
                [ke.Ma, ke.Base] = ki.interpolate(fxi);
                
                fsvr = ke.evaluate(x);
                
                h = pm.nextPlot('xx',sprintf('Interpolation test with kernel %s, newton=%d',class(ke.Kernel),lu));
                hold(h,'on');
                for k=1:size(fx,1)
                    plot(h,x,fx(k,:),'r');
                    % Plot approximated function
                    plot(h,x,fsvr(k,:),'b--');
                    %skipped = setdiff(1:length(x),svidx);
                    plot(h,xi,fxi(k,:),'.r');
                end
            end
        end
        
        function test_KernelInterpolationError
            realprec = 20;
            range = [-4 4];
            
            pm = PlotManager(false,2,3);
            pm.LeaveOpen = true;
            
            ke = kernels.KernelExpansion;
            ke.Kernel = kernels.GaussKernel(1);
            testfun = @(x)sinc(x) * 2 .*sin(3*x);
            ke.Kernel = kernels.InvMultiquadrics(-1,1);
            
            ki = general.interpolation.KernelInterpol;
            
            hsteps = [.5:.25:2 2:2:10];
            err = zeros(size(hsteps));
            for n = 1:length(hsteps)
                h = 1/hsteps(n);
                xi = range(1):h:range(2);
                xi = xi + rand(size(xi))*max(xi)/100+.1;
                fxi = testfun(xi);
                
                x = range(1):h/realprec:range(2);
                fx = testfun(x);
                
                ke.Centers.xi = xi;
                ki.init(ke);
                [ke.Ma, ke.Base, used] = ki.interpolate(fxi);
                ke.Centers.xi = xi(:,used);
                fint = ke.evaluate(x);
                
                ax = pm.nextPlot(['step' n],sprintf('h=%g',h),'x','fx');
                plot(ax,xi,fxi,'r.',x,fx,'r',x,fint,'b');
                err(n) = max(abs(fx-fint));
            end
            pm = PlotManager; pm.LeaveOpen = true;
            ax = pm.nextPlot('errors','Error development over h sizes','1/h','errors');
            h = 1./hsteps;
            semilogy(ax,hsteps,exp(log(h)./h),'b',hsteps,exp(-1./h),'b--',hsteps,err,'r');
            legend(ax,'exp(log(h)/h)','exp(-1/h)','||f-sfx||_\infty');
            pm.done;
        end
    end
end

