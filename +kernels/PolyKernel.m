classdef PolyKernel < kernels.BaseKernel
    %POLYKERNEL Basic polynomial kernel
    %
    % @author Daniel Wirtz @date 12.03.2010
    
    properties
        % The degree of the polynomial kernel.
        %
        % Defaults to 2.
        Degree = 2;
    end
    
    methods
        
        function this = PolyKernel(deg)
            if nargin == 1
                % TODO validity checks
                this.Degree = deg;
            end
            this.RotationInvariant = true;
        end
        
        function c = getGlobalLipschitz(this)
            error('Not implemented yet!');
        end
        
        function K = evaluate(this, x, y)
            if nargin == 2
                y=x;
            end
            K = x'*y.^this.Degree;
        end
    end
    
end
