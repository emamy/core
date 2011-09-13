classdef KernelInterpol < KerMorObject & approx.algorithms.IKernelCoeffComp
    % Provides kernel interpolation.
    %
    % The basic interpolation form is 
    % `` f(x) = \sum\limits_{i=1}^N \alpha_i \Phi(x,x_i)``
    % Interpolation finds coefficients such that 
    % `fx_i = RHS(x_i)` for `i=1\ldotsN`.
    %
    % There is also a zero-function threshold `10*eps`. If all
    % `fx_i-\beta` values are below that a constant function is
    % assumed.
    %
    % @author Daniel Wirtz @date 01.04.2010
    %
    % @change{0,5,dw,2011-09-12} Set the UseLU flag to true per default.
    % Using AKernelMatrix instances now, along with flags of whether to
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
        K;
    end
    
    properties(Access=private)
        fK;
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

%             if all(abs(fxi) < 10*eps)
%                 a = zeros(size(fxi))';
%                 if KerMor.App.Verbose > 3
%                     fprintf('KernelInterpol note: All fxi values < 10eps, assuming zero coefficients!\n');
%                 end
%             else
            if isa(this.fK, 'data.MemoryKernelMatrix')
                if this.fK.UseLU
                    a = this.fK.U\(this.fK.L\fxi');
                    return;
                elseif this.fK.BuildInverse
                    a = this.fK.Kinv * fxi';
                    return;
                end
            end
            a = this.fK\fxi';
%             end
        end
        
        function set.K(this, value)
            % Sets the kernel matrix property and computes the LU
            % decomposition if the UseLU flag is set to true.
            if ~isa(value, 'data.AKernelMatrix')
                error('K must be a data.AKernelMatrix');
            end
            this.fK = value;
        end
        
        function K = get.K(this)
            K = this.fK;
        end
        
        %% approx.algorithms.IKernelCoeffComp interface members
        function init(this, K)
            this.K = K;
        end
        
        function [ai, svidx] = computeKernelCoefficients(this, yi, dummy)%#ok
            
            % Transform to row vector
            ai = this.interpolate(yi)';
            svidx = 1:size(ai,2);
        end
    end
    
    methods(Static)
        
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
            
            %x = -5:.05:5;
            %fx = sinc(x);
            
            x = 0:.05:10;
            fx = x.^2.*sin(x);
            
            n = length(x);
            samp = 1:15:n;
            
            k = kernels.GaussKernel(1);
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

