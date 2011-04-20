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
        kappa = 1;
        
        % The constant inside tanh
        nu = -1;
    end
    
    methods
        
        function this = SigmoidKernel(kappa, nu)
            if nargin > 0
                this.kappa = kappa;
                if nargin > 1
                    this.nu = nu;
                end
            end
            this.RotationInvariant = true;
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
            if nargin == 2
                y = x;
            end
            K = tanh( this.kappa * x'*y + this.nu);
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

