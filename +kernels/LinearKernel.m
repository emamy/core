classdef LinearKernel < kernels.BaseKernel
    %LINEARKERNEL The simple scalar-product kernel.
    %   Detailed explanation goes here
    
    methods
        
        function c = getGlobalLipschitz(this)%#ok
            c = 1;
        end
        
        function Nabla = getNabla(this, x, y)%#ok
            % Partial derivatives of scalar product is simply the second argument vector.
            Nabla = y;
        end
        
        function K = evaluate(this, x, y)%#ok
            if nargin == 2
                y = x;
            end
            K = x'*y;
        end
    end
    
end

