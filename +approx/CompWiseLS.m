classdef CompWiseLS < approx.BaseCompWiseKernelApprox
    %COMPWISELS Summary of this class goes here
    %   Detailed explanation goes here
    
    %| @docupdate
    
    properties
        % The `\lambda` least squares constant. Forwarded to the single
        % dimension kernel least squares algorithm.
        %
        % See also: general.regression.KernelLS
        Lambda=1;
    end
    
    properties(Access=private)
        % The used single-dimension kernel least-squares algorithm.
        %
        % @type general.regression.KernelLS
        %
        % See also: general.regression.KernelLS
        LS;
    end
    
    methods(Sealed)
        function target = clone(this)
            % Makes a copy of this instance.
            %
            % See also: ICloneable
            
            target = approx.CompWiseLS;
            
            % Copy local properties
            target.Lambda = this.Lambda;
            % Clone KernelLS instance
            ls = general.regression.KernelLS;
            ls.lambda = this.LS.lambda;
            ls.CGMaxIt = this.LS.CGMaxIt;
            ls.CGTol = this.LS.CGTol;
            ls.MaxStraightInvDim = this.LS.MaxStraightInvDim;
            target.LS = ls;
            
            % Call superclass clone
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
        end
    end
    
    methods(Access=protected, Sealed)
        
        function prepareApproximationGeneration(this, K)
            % @copydoc approx.BaseCompWiseKernelApprox
            
            this.LS = general.regression.KernelLS;
            % Assign Kernel matrix
            this.LS.K = K;
            this.LS.lambda = this.Lambda;
        end
        
        function [ai, b, svidx] = calcComponentApproximation(this, fxi)
            ai = this.LS.regress(fxi);
            b = 0;
            svidx = [];
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
            m.Data.ApproxTrainData = [zeros(2,size(x,2)); m.Times; x];
            m.Data.ApproxfValues = fxi;
            a = approx.CompWiseLS;
            a.Lambda = .1;
            a.TimeKernel = kernels.LinearKernel;
            a.SystemKernel = kernels.GaussKernel(5);
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
            m.Data.ApproxTrainData = [zeros(2,size(x,2)); m.Times; x];
            m.Data.ApproxfValues = fxi;
            a = approx.CompWiseLS;
            a.Lambda = 1;
            
            a.TimeKernel = kernels.GaussKernel;
            %a.TimeKernel = kernels.LinearKernel;
            
            c = kernels.CombinationKernel;
            c.addKernel(kernels.LinearKernel,10);
            %c.addKernel(kernels.PolyKernel(2));
            c.addKernel(kernels.PolyKernel(3));
            %c.addKernel(kernels.PolyKernel(4));
            
            a.SystemKernel = c;
            %a.SystemKernel = kernels.GaussKernel(4);
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

