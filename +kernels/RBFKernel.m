classdef RBFKernel < kernels.BaseKernel
    %RBFKERNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Gamma = 1;
    end
    
    methods
        function this = RBFKernel(Gamma)
            if nargin == 1
                this.Gamma = Gamma;
            end
            this.RotationInvariant = true;
        end
        
        function K = evaluate(this, x, y)
            % Original source from Bernard Haasdonk / KerMet
            n1sq = sum(x.^2,1);
            n1 = size(x,2);
            
            if nargin == 2;
                n2sq = n1sq;
                n2 = n1;
                y = x;
            else
                n2sq = sum(y.^2,1);
                n2 = size(y,2);
            end;
            K = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*x'*y;
            K(K<0) = 0;
            K = exp(-this.Gamma*K);
        end
    end
    
end

