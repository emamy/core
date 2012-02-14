classdef PolyKernel < kernels.BaseKernel
    %POLYKERNEL Basic polynomial kernel
    %
    % @author Daniel Wirtz @date 12.03.2010
    
    properties(SetObservable)
        % The degree of the polynomial kernel.
        %
        % @propclass{critical} Greatly influences the kernels behaviour.
        %
        % @default 2
        Degree = 2;
    end
    
    methods
        
        function this = PolyKernel(deg)
            if nargin == 1
                % TODO validity checks
                this.Degree = deg;
            end
        end
        
        function c = getGlobalLipschitz(this)%#ok
            % @todo implement
            error('Not implemented yet!');
        end
        
        function Nabla = getNabla(this, x, y)
            % Partial derivatives of scalar product is simply the second argument vector.
            Nabla = this.Degree*bsxfun(@times, y, x'*y.^(this.Degree-1));
        end
        
        function K = evaluate(this, x, y)
            if nargin == 2 || isempty(y)
                y=x;
            end
            K = x'*y.^this.Degree;
        end
    end
    
end

