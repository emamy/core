classdef LinearKernel < kernels.BaseKernel
    %LINEARKERNEL The simple scalar-product kernel.
    %   Detailed explanation goes here
    
    methods
        
        function this = LinearKernel
            this.RotationInvariant = true;
        end
        
        
        function K = evaluate(this, x, y)
            if nargin == 2
                y = x;
            end
            K = x'*y;
        end
    end
    
end

