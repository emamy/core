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
                % \todo validity checks
                this.Degree = deg;
            end
            this.IsScProd = true;
        end
        
        function copy = clone(this)
            copy = clone@kernels.BaseKernel(this, kernels.PolyKernel);
            copy.Degree = this.Degree;
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
            if ~isempty(this.fP)
                x = x(this.fP,:);
            end
            if nargin == 2 || isempty(y)
                y = x;
            else
                if ~isempty(this.fP)
                    y = y(this.fP,:);
                end
            end
            K = (x'*(this.fG*y)).^this.Degree;
%             sx = this.fG*x;
%             n1sq = sum(x.*sx,1);
%             n1 = size(x,2);
%             if nargin == 2 || isempty(y)
%                 y = x;
%                 n2sq = n1sq;
%                 n2 = n1;
%             else
%                 if ~isempty(this.fP)
%                     y = y(this.fP,:);
%                 end
%                 n2sq = sum(y.*(this.fG*y),1);
%                 n2 = size(y,2);
%             end
%             no = sqrt(repmat(n1sq',1,n2) .* repmat(n2sq,n1,1));
%             a = (sx'*y);
%             b = a ./ no;
%             K = b.^this.Degree;
%             %K = ((sx'*y) ./ no).^this.Degree;
        end
    end
    
end

