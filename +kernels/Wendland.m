classdef Wendland < kernels.ARBFKernel
    % Wendland: Implementation of various Wendland kernels
    %
    % Further background e.g. in §9.4 of @cite{W05}, especially Corollary 9.14.
    %
    % @author Daniel Wirtz @date 2013-01-16
    %
    % @new{0,7,dw,2013-01-16} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Dependent)
        % The considered space dimension
        d;
        
        % Desired smoothness `2k`
        k;
    end
    
    properties(Access=private)
        fd = 1;
        fk = 0;
        expo;
        co;
        polystr;
        polyfun;
    end
    
    methods
        function this = Wendland
            this.polystr = {...
                '@(r)(l+1)*r+1' ...
                '@(r)(l^2+4*l+3)*r.^2/3+(l+2)*r+1' ...
                '@(r)((l^3+9*l^2+23*l+15)*r.^3 + (6*l^2+36*l+45)*r.^2)/15+(l+3)*r+1'};
            this.updateCoeffs;
        end
        
        function Kxy = evaluate(this, x, y)
            r = sqrt(this.getSqDiffNorm(x, y))/this.Gamma;
            rp = max(1-r, 0);
            %rp = (1 - r).*(r <= 1);
            
            p = 1;
            if (this.fk > 0)
                p = this.polyfun(r);
            end
            Kxy=(rp.^this.expo).*p;
        end
        
        function Kxy = evaluateScalar(this, r)
            error('not implemented');
            rp = max(1-r,0);
            %rp = (1-r).*(r <= 1);
            
            p = 1;
            if (this.fk > 0)
                p = this.polyfun(r);
            end
            Kxy=(rp.^this.expo).*p;
        end 
        
        function Nabla = getNabla(this, x, y)
            error('not implemented');
        end
        
        function dc = getDefaultConfig(this)
            dc = kernels.config.WendlandConfig('S',this.k,'Dim',this.d);
        end
        
        % Returns the global lipschitz constant of this kernel.
        %
        % Exprimental state as not implemented & checked for all kernels.
        function c = getGlobalLipschitz(this)
            error('not implemented');
        end
        
        function copy = clone(this)
            copy = clone@kernels.ARBFKernel(this, kernels.Wendland);
            % Triggers update of coeffs
            copy.d = this.fd;
            copy.k = this.fk;
        end
    end
    
    methods(Access=private)
        function updateCoeffs(this)
            % calculates coefficients and exponent for polynomial
            % part p_{d,k}(r) of the Wendland kernel function.
            l = floor(this.fd/2) + this.fk + 1;
            this.expo = l + this.fk;
            this.polyfun = [];
            if (this.fk > 0)
                this.polyfun = eval(this.polystr{this.fk});
            end
        end
    end
    
    %% Getter & Setter
    methods
        function set.k(this, value)
            if value < 0 || value > 3
                error('Only k values in [0,3] are allowed');
            elseif round(value) ~= value
                error('Only the k values 0,1,2,3 are allowed');
            end
            this.fk = value;
            this.updateCoeffs;
        end
        
        function set.d(this, value)
            if round(value) ~= value
                error('The dimension must be an integer');
            end
            this.fd = value;
            this.updateCoeffs;
        end
        
        function v = get.d(this)
            v = this.fd;
        end
        
        function v = get.k(this)
            v = this.fk;
        end
    end
    
    methods(Static)
        function res = test_WendlandKernel(pm)
            if nargin < 1
                pm = PlotManager(false,3,3);
                pm.LeaveOpen = true;
            end
            c = 0;
            x = (-1.2:.01:1.2)+c;
            [X,Y] = meshgrid(x);
            x2 = [X(:)'; Y(:)'];
            k = kernels.Wendland;
            kexp = kernels.KernelExpansion;
            kexp.Kernel = k;
            kexp.Centers.xi = c;
            kexp.Ma = 1;
            conf = general.Utils.createCombinations([1 2 3 4 5],[0 1 2 3]);
            for n = 1:length(conf)
                k.d = conf(1,n);
                k.k = conf(2,n);
                tag = sprintf('w_1d_d%d_k%d',k.d,k.k);
                h = pm.nextPlot(tag,sprintf('Wendland kernel with d=%d,k=%d on 1D data',k.d,k.k));
                plot(h,x,kexp.evaluate(x));
            end
            kexp.Centers.xi = [c; c];
            for n = 1:length(conf)
                k.d = conf(1,n);
                k.k = conf(2,n);
                tag = sprintf('w_2d_d%d_k%d',k.d,k.k);
                h = pm.nextPlot(tag,sprintf('Wendland kernel with d=%d,k=%d on 2D data',k.d,k.k));
                surf(h,X,Y,reshape(kexp.evaluate(x2),length(x),[]),'EdgeColor','none');
            end
            if nargin < 1
                pm.done;
            end
            res = true;
        end
    end
end