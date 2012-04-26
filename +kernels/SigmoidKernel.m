classdef SigmoidKernel < kernels.BaseKernel
    %SIGMOIDKERNEL The sigmoid kernel.
    %
    % This kernel is defined via `K(x,y) = \tanh(\kappa<x,y> + \nu)`
    %
    % See also: BaseKernel LinearKernel PolyKernel GaussKernel
    % CombinationKernel
    %
    % @author Daniel Wirtz 23.03.2010
    
    properties
        % The factor for the scalar product inside tanh
        %
        % @propclass{critical} Greatly influences the kernels behaviour.
        kappa = 1;
        
        % The constant inside tanh
        %
        % @propclass{critical} Greatly influences the kernels behaviour.
        nu = -1;
    end
    
    methods
        
        function this = SigmoidKernel(kappa, nu)
            this.registerProps('kappa','nu');
            this.IsScProd = true;
            if nargin > 0
                this.kappa = kappa;
                if nargin > 1
                    this.nu = nu;
                end
            end
        end
        
        function copy = clone(this)
            copy = clone@kernels.BaseKernel(this, kernels.SigmoidKernel);
            copy.kappa = this.kappa;
            copy.nu = this.nu;
        end
        
        function c = getGlobalLipschitz(this)%#ok
            % @todo implement
            error('Not implemented yet!');
        end
        
        function Nabla = getNabla(this, x, y)%#ok
            % @todo implement
            error('Not implemented yet!');
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
            K = tanh( this.kappa * x'*(this.G*y) + this.nu);
        end
        
        function set.kappa(this, value)
            if ~isposrealscalar(value)
                error('kappa must be a real scalar and positive');
            end
            this.kappa = value;
        end
        
        function set.nu(this, value)
            if ~isposrealscalar(-value)
                error('kappa must be a real scalar and negative');
            end
            this.nu = value;
        end
    end
    
end

