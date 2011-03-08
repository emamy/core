classdef CompWiseSVR < approx.BaseCompWiseKernelApprox
    % Performs component-wise support vector regression.
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    properties
        % The SVR instance to use. Must be a descendend of
        % general.regression.BaseScalarSVR
        %
        % Default: general.regression.ScalarNuSVR
        ScalarSVR;
    end

    methods
        function this = CompWiseSVR
            this.ScalarSVR = general.regression.ScalarNuSVR;
            this.ScalarSVR.nu = .5;
        end
    end
    
    methods(Sealed)
        function target = clone(this)
            % Makes a copy of this instance.
            %
            % See also: ICloneable
            
            target = approx.CompWiseSVR;
            % Call superclass clone
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
            % Copy local properties
            target.ScalarSVR = this.ScalarSVR.clone;
        end
    end
    
    methods(Access=protected, Sealed)
        
        function prepareApproximationGeneration(this, K)
            this.ScalarSVR.K = K;
        end
        
        function [ai,b,svidx] = calcComponentApproximation(this, fxi)
            [ai,b,svidx] = this.ScalarSVR.regress(fxi);
        end

    end
    
    methods(Static)
        
        function res = test_Cloning
            a = approx.CompWiseSVR;
            a2 = a.clone;
            res = ~a.eq(a2);
            clear a a2;
        end
        
        function test_CompWiseSVR2D2D
            dim = 3;
            [X,Y] = meshgrid(-1:1/dim:1,-1:1/dim:1);
            
            fx1 = sin(pi*X) + Y;
            fx2 = .5*exp(abs(X-Y));
            fx3 = sin(X.*Y);
            fx4 = X.^2+5*cos(Y);
            
            xi(1,:) = X(:);
            xi(2,:) = Y(:);
            fxi(1,:) = fx1(:);
            fxi(2,:) = fx2(:);
            fxi(3,:) = fx3(:);
            fxi(4,:) = fx4(:);
            
            
            m = models.BaseFullModel;
            m.T = size(xi,2)-1;
            m.dt = 1;
            m.Data.ApproxTrainData = [zeros(2,size(xi,2)); m.scaledTimes; xi];
            m.Data.ApproxfValues = fxi;
            a = approx.CompWiseSVR;
            %a.ScalarSVR = general.regression.ScalarEpsSVR;
            %a.ScalarSVR.eps = .1;
            a.ScalarSVR.C = 100000;
            a.TimeKernel = kernels.GaussKernel(4);
            a.SystemKernel = kernels.GaussKernel(4);
            m.Approx = a;
            m.Approx.approximateCoreFun(m);
            fxiap = m.Approx.evaluate(xi,m.scaledTimes,[]);
            
            figure(1)
            subplot(2,2,1);
            surfl(X,Y,fx1);
            hold on;
            mesh(X,Y,reshape(fxiap(1,:),size(X,1),[]));
            
            subplot(2,2,2);
            surfl(X,Y,fx2);
            hold on;
            mesh(X,Y,reshape(fxiap(2,:),size(X,1),[]));
            
            subplot(2,2,3);
            surfl(X,Y,fx3);
            hold on;
            mesh(X,Y,reshape(fxiap(3,:),size(X,1),[]));
            
            subplot(2,2,4);
            surfl(X,Y,fx4);
            hold on;
            mesh(X,Y,reshape(fxiap(4,:),size(X,1),[]));
            
        end
        
        function test_CompWiseSVR1D4D
            dim = 10;
            X = -1:1/dim:1;
            
            fx1 = sin(pi*X);
            fx2 = .5*exp(X);
            fx3 = sin(5*X.^2);
            fx4 = X.^2+5*cos(X);
            
            fxi(1,:) = fx1;
            fxi(2,:) = fx2;
            fxi(3,:) = fx3;
            fxi(4,:) = fx4;
            
            
            m = models.BaseFullModel;
            m.T = size(X,2)-1;
            m.dt = 1;
            m.Data.ApproxTrainData = [zeros(2,size(X,2)); m.scaledTimes;X];
            m.Data.ApproxfValues = fxi;
            a = approx.CompWiseSVR;
            a.ScalarSVR = general.regression.ScalarEpsSVR;
            a.ScalarSVR.eps = .1;
            a.ScalarSVR.C = 100000;
            c = kernels.CombinationKernel;
            c.addKernel(kernels.GaussKernel(1));
            c.addKernel(kernels.GaussKernel(6));
            c.addKernel(kernels.LinearKernel);
            c.CombOp = @(x,y)x.*y;
            a.SystemKernel = c;
            a.TimeKernel = kernels.GaussKernel(7);
            m.Approx = a;
            m.Approx.approximateCoreFun(m);
            fxiap = m.Approx.evaluate(X,m.scaledTimes,[]);
            
            figure(2);
            subplot(2,2,1);
            plot(X,fx1,'r',X,fxiap(1,:),'b-');
            
            subplot(2,2,2);
            plot(X,fx2,'r',X,fxiap(2,:),'b-');
            
            subplot(2,2,3);
            plot(X,fx3,'r',X,fxiap(3,:),'b-');
            
            subplot(2,2,4);
            plot(X,fx4,'r',X,fxiap(4,:),'b-');
            
        end
    end
    
end