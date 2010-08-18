classdef KernelInterpol
    %KERNELINTERPOL This is the class summary.
    % Some detailed information!
    % Yeay a texty text. Again, a texty text.
    %
    % @ingroup general_interpolation
    %
    % @author Daniel Wirtz
    
    properties
        K;
    end
    
    methods
        function a = interpolate(this, fxi)
            % Interpolates the 
            a = this.K\(fxi');
        end
    end
    
    methods(Static)
        function test_KernelInterpolation()
            % Performs a test of this class
            
            x = -5:.05:5;
            fx = sinc(x);
            
            n = length(x);
            samp = 1:15:n;
            
            k = kernels.GaussKernel(1);
            internal_test(x,fx,samp,k,10);
            
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

