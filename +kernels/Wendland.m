classdef Wendland < kernels.ARBFKernel
    % Wendland: Implementation of Wendland kernel `phi_{d,k}(r)=p_{d,k}(r)*(1-r)^e`
    %
    % Further background e.g. in §9.4 of @cite{W05}.
    %
    % Original code from R. Schaback etc
    % We use only those which are pos. def. in dimension at most 3.
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
    
    properties%(Access=private)
        fd = 1;
        fk = 0;
        expo;
        co;
        polystr;
        polyfun;
    end
    
    methods
        function this = Wendland
            this.updateCoeffs;
            this.polystr = {...
                '@(r)(l+1)*r+1' ...
                '@(r)(l^2+4*l+1)*r.^2+(3*l+6)*r+1' ...
                '@(r)(l^3+9*l^2+25*l+15)*r.^3 + (6*l^2+36*l+45)*r.^2+(15*l+45)*r+15'};
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
        
        function phir = evaluateScalar(this, r)
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
        
        % Returns the global lipschitz constant of this kernel.
        %
        % Exprimental state as not implemented & checked for all kernels.
        function c = getGlobalLipschitz(this)
            error('not implemented');
        end
    end
    
    methods(Access=private)
        function updateCoeffs(this)
            % calculates coefficients and exponent for polynomial
            % part p_{d,k}(r) of the Wendland kernel function.
            l = floor(this.fd/2) + this.fk + 1;
            this.expo = l + (this.fk > 0);
            this.polyfun = [];
            if (this.fk > 0)
                this.polyfun = eval(this.polystr{this.fk});
            end
%             expon = floor(this.fd/2)+this.fk+1;
%             coeff = zeros(this.fk+1,1);
%             coeff(1,1)=1;
%             for n=0:this.fk-1
%                 coeff(n+2,1)=coeff(n+1,1)/(n+expon+2);
%                 for j=n+1:-1:2
%                     coeff(j,1)=(j*coeff(j+1,1)+coeff(j-1,1))/(expon+j);
%                 end
%                 expon=expon+1;
%                 coeff(1,1)=coeff(2,1)/expon;
%             end
%             this.expo = expon;
%             this.co = coeff;
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
end