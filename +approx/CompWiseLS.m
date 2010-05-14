classdef CompWiseLS < approx.BaseKernelApprox
    %COMPWISELS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda=1;
    end
    
    properties(Access=private)
        adata;
        xi;
    end
    
    methods
        function this = CompWiseLS
            this.CustomProjection = true;
        end
    end
    
    methods(Access=protected)
        
        function fx = evaluate_approximation(this, x)
            % Evaluates the approximated function at point x
            K = this.evaluateKernel(x,this.xi);
            fdims = size(this.adata,1);
            fx = zeros(fdims,size(x,2));
            for fdim = 1:fdims
                fx(fdim,:) = K * this.adata(fdim,:)';
            end
        end
        
        function gen_approximation_data(this, xi, fxi)
            % Computes the approximation according to the concrete
            % approximation strategy. 
            
            ls = general.regression.KernelLS;
            ls.lambda = this.lambda;
            ls.K = this.evaluateKernel(xi);
            
            % Make hermetian (no rounding errors)
            %ls.K = .5*(ls.K'+ls.K);
            
            wh = waitbar(0,'Initializing component-wise kernel LS');
            try
                fdims = size(fxi,1);
                this.adata = zeros(fdims, size(xi,2));
                for fdim = 1:fdims
                    waitbar(fdim/fdims,wh,sprintf('Performing LS for dimension %d/%d ... %2.0f %%',fdim,fdims,(fdim/fdims)*100));
                    fx = fxi(fdim,:);
                    this.adata(fdim,:) = ls.regress(fx);
                end
                close(wh);
            catch ME
                close(wh);
                rethrow(ME)
            end
            this.xi = xi;
        end
        
        function copy = customProject(this, V)
            copy = this.clone;
            
            copy.adata = V'*this.adata;
            
            if this.RotationInvariantKernel
                % Extract system part and project into V space
                [x,t,mu] = this.splitTripleVect(this.xi);
                x = V' * x;
                copy.xi = this.compileTripleVect(x,t,mu);
            end
        end
        
         function target = clone(this)
            % Makes a copy of this instance.
            %
            % See also: ICloneable
            
            target = approx.CompWiseLS;
            % Call superclass clone
            target = clone@approx.BaseKernelApprox(this, target);
            % Copy local properties
            target.lambda = this.lambda;
            
            target.adata = this.adata;
            target.xi = this.xi;
        end
    end
    
    methods(Static)
        function test_CompWiseLS2D
            dim = 6;
            [X,Y] = meshgrid(-1:1/dim:1,-1:1/dim:1);
            %x = -1:1/dim:1;
            
            fx1 = sin(pi*X) + Y;
            fx2 = .5*exp(abs(X-Y));
            
            x(1,:) = X(:);
            x(2,:) = Y(:);
            fxi(1,:) = fx1(:);
            fxi(2,:) = fx2(:);
            
            m = models.BaseFullModel;
            m.T = size(x,2)-1;
            m.dt = 1;
            m.Data.Snapshots = x;
            m.Data.fValues = fxi;
            a = approx.CompWiseLS;
            a.lambda = .1;
            a.TimeKernel = kernels.LinearKernel;
            a.SystemKernel = kernels.RBFKernel(5);
            m.Approx = a;
            m.Approx.approximateCoreFun(m);
            
            [X2,Y2] = meshgrid(-1:.5/dim:1,-1:.5/dim:1);
            x2(1,:) = X2(:);
            x2(2,:) = Y2(:);
            % Fit t (insert "mid-times")
            m.dt = m.T/(size(x2,2)-1);
            fxiap = m.Approx.evaluate(x2,m.Times,[]);
            
            figure(1)
            subplot(1,2,1);
            surfl(X,Y,fx1);
            hold on;
            mesh(X2,Y2,reshape(fxiap(1,:),size(X2,1),[]));
            
            subplot(1,2,2);
            surfl(X,Y,fx2);
            hold on;
            mesh(X2,Y2,reshape(fxiap(2,:),size(X2,1),[]));
            
        end
        
        function test_CompWiseLS1D
            dim = 60;
            x = -2:1/dim:2;
            
            fxi(1,:) = x+sin(pi*x);
            %fxi(2,:) = .5*exp(abs(x));
            fxi(2,:) = cos(.5*exp(x));
            
            m = models.BaseFullModel;
            m.T = length(x)-1;
            m.dt = 1;
            m.Data.Snapshots = x;
            m.Data.fValues = fxi;
            a = approx.CompWiseLS;
            a.lambda = 1;
            
            a.TimeKernel = kernels.RBFKernel;
            %a.TimeKernel = kernels.LinearKernel;
            
            c = kernels.CombinationKernel;
            c.addKernel(kernels.LinearKernel,10);
            %c.addKernel(kernels.PolyKernel(2));
            c.addKernel(kernels.PolyKernel(3));
            %c.addKernel(kernels.PolyKernel(4));
            
            a.SystemKernel = c;
            %a.SystemKernel = kernels.RBFKernel(4);
            %a.SystemKernel = kernels.PolyKernel(3);
            
            m.Approx = a;
            m.Approx.approximateCoreFun(m);
            
            x2 = -2:.5/dim:2;
            % Fit t (insert "mid-times")
            m.dt = m.T/(length(x2)-1);
            fxiap = m.Approx.evaluate(x2,m.Times,[]);
            
            figure(1)
            subplot(1,2,1);
            plot(x,fxi(1,:),'r',x2,fxiap(1,:),'b');
            legend('full','approx');
            axis tight;
            
            subplot(1,2,2);
            plot(x,fxi(2,:),'r',x2,fxiap(2,:),'b');
            legend('full','approx');
            axis tight;
            
        end
    end
    
end

