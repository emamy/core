classdef LinearKernel < kernels.BaseKernel
    %LINEARKERNEL The simple scalar-product kernel.
    %   Detailed explanation goes here
    
    methods
        function this = LinearKernel
            this.IsScProd = true;
        end
        
        function c = getGlobalLipschitz(this)%#ok
            c = 1;
        end
        
        function Nabla = getNabla(this, x, y)%#ok
            % Partial derivatives of scalar product is simply the second argument vector.
            Nabla = y;
        end
        
        function K = evaluate(this, x, y)
            if ~isempty(this.P)
                x = x(this.P,:);
            end
            if nargin == 2 || isempty(y)
                y = x;
            else
                if ~isempty(this.P)
                    y = y(this.P,:);
                end
            end
            K = x'*(this.G*y);
        end
        
        function copy = clone(this)
            copy = clone@kernels.BaseKernel(this, kernels.LinearKernel);
        end
    end
    
end

