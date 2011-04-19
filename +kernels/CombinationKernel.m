classdef CombinationKernel < kernels.BaseKernel
    %SUMKERNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The combination operation for each kernel.
        %
        % Has to be a function handle taking two arguments:
        % first: The current accumulation result
        % second: The last kernel's evaluation
        CombOp = @(a,b)a+b;
        
        %CombOp = @(a,b)a.*b;
        %CombOp = @(a,b)a.^b;
    end
    
    properties(Access=private)
        Kernels = {};
        Weights = [];
    end
    
    methods
        
        function this = CombinationKernel
            % Invariance true as long as no kernel is assigned :-)
            this.RotationInvariant = true;
        end
        
        function K = evaluate(this, x, varargin)
            if isempty(this.Kernels)
                error('No Kernels here to combine. Use addKernel for setup.');
            end
            K = this.Weights(1)*this.Kernels{1}.evaluate(x,varargin{:});
            for n = 2:length(this.Kernels)
                K = this.CombOp(K, this.Weights(n)*this.Kernels{n}.evaluate(x,varargin{:}));
            end
        end
        
        function c = getGlobalLipschitz(this)%#ok
            % @todo implement
            error('Not implemented yet');
        end
        
        function Nabla = getNabla(this, x, y)%#ok
            % @todo implement
            error('Not implemented yet');
        end
        
        function this = addKernel(this, kernel, weight)
            if nargin == 2
                weight = 1;
            end
            if ~isa(kernel,'kernels.BaseKernel')
                error('kernel param must be a kernels.BaseKernel!');
            end
            this.Kernels{end+1} = kernel;
            this.Weights(end+1) = weight;
            this.RotationInvariant = this.RotationInvariant && kernel.RotationInvariant;
        end
    end
    
    methods(Static)
        function test_CombinationKernels
            k1 = kernels.CombinationKernel;
            k1.addKernel(kernels.GaussKernel(.7));
            k1.addKernel(kernels.GaussKernel(.5));
            %k1.addKernel(kernels.PolyKernel(3),3);
            k1.addKernel(kernels.PolyKernel(2),2);
            k1.CombOp = @(a,b)a.*b;
            k2 = kernels.CombinationKernel;
            k2.addKernel(kernels.LinearKernel,4);
            k2.addKernel(kernels.PolyKernel(2));
            k2.CombOp = @(a,b)a+b;
            k = kernels.CombinationKernel;
            k.addKernel(k1);
            k.addKernel(k2);
            
            x = -1:.05:1;
            fx = sinc(x);
            
            svr = general.regression.ScalarEpsSVR;
            svr.QPSolver = solvers.qp.qpOASES;
            svr.eps = .3;
            svr.C = 2;
            svr.K = k.evaluate(x,x);
            
            figure(1);
            plot(x,fx,'r',x,[fx-svr.eps; fx+svr.eps],'r--');
            
            [ai,b,svidx] = svr.regress(fx);
            sv = x(svidx);
            svfun = @(x)ai'*k.evaluate(sv,x) + b;
            
            fsvr = svfun(x);
            
            hold on;
            
            % Plot approximated function
            plot(x,fsvr,'b',x,[fsvr-svr.eps; fsvr+svr.eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');
            
            hold off;
        end
    end
    
end
