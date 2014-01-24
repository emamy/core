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
    
    methods(Static)
        function res = test_PolyKernel(pm)
            if nargin < 1
                pm = PlotManager(false,3,3);
                pm.LeaveOpen = true;
            end
            c = 1;
            x = (-1.2:.01:1.2)+c;
            
            k = kernels.PolyKernel;
            kexp = kernels.KernelExpansion;
            kexp.Kernel = k;
            kexp.Centers.xi = c;
            kexp.Ma = 1;
            conf = [.1 .4 .8 1 2 3 4];
            for n = 1:length(conf)
                k.Degree = conf(n);
                tag = sprintf('poly_1d_deg%g',k.Degree);
                h = pm.nextPlot(tag,sprintf('Polynomial kernel with deg=%g on 1D data',k.Degree));
                z = kexp.evaluate(x);
                plot(h,x,[real(z); imag(z)]);
            end
            kexp.Centers.xi = [c; 1.6*c];
            [X,Y] = meshgrid(x);
            x2 = [X(:)'; Y(:)'];
            for n = 1:length(conf)
                k.Degree = conf(n);
                tag = sprintf('poly_2d_deg%g',k.Degree);
                h = pm.nextPlot(tag,sprintf('Polynomial kernel with deg=%g on 2D data',k.Degree));
                Z = real(reshape(kexp.evaluate(x2),length(x),[]));
                surf(h,X,Y,Z,'EdgeColor','none');
            end
            if nargin < 1
                pm.done;
            end
            res = true;
        end
    end
    
end

