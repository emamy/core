classdef NoKernel < kernels.ARBFKernel
    % Neutral Kernel which has no effect.
    %
    % A call to evaluate just returns 1.
    
    methods
        
        function this = NoKernel
            this = this@kernels.ARBFKernel;
            this.G = [];
            this.P = [];
        end
        
        function copy = clone(this)
            copy = clone@kernels.ARBFKernel(this, kernels.NoKernel);
        end
        
        function K = evaluate(this, x, y)%#ok
            K = 1;
        end
        
        function K = evaluateScalar(this, x)%#ok
            K = 1;
        end
        
        function Nablax = getNabla(this, x, y)%#ok
            % Return zero as this kernel is a constant one.
            warning('KerMor:uninvestigated','Using this function may cause unpredictable behaviour.');
            Nablax = 0;
        end
        
        function c = getGlobalLipschitz(this)%#ok
            c = 1;
        end
    end
end

