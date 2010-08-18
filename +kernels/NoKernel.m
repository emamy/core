classdef NoKernel < kernels.BaseKernel
    % Neutral Kernel which has no effect.
    %
    % A call to evaluate just returns 1.
    
    methods
        
        function this = NoKernel
            this.RotationInvariant = true;
        end
        
        function K = evaluate(this, x, y)%#ok
            K = 1;
        end
        
        function c = getGlobalLipschitz(this)%#ok
            c = 1;
        end
    
    end
    
end

