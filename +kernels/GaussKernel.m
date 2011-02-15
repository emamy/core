classdef GaussKernel < kernels.BellFunction
    %RBFKERNEL Summary of this class goes here
    %   Radial Basis Function Kernel.
    %   Uses the notation
    %   ``\Phi(x,y) = e^{\frac{||x-y||^2}{\gamma}},``
    %   so be careful with the `\gamma` constant.
    
    properties
        Gamma = 1;
    end
    
    methods
        function this = GaussKernel(Gamma)
            if nargin == 1
                this.Gamma = Gamma;
            end
        end
        
        function K = evaluate(this, x, y)
            n1sq = sum(x.^2,1);
            n1 = size(x,2);
            
            if nargin == 2;
                n2sq = n1sq;
                n2 = n1;
                y = x;
            else
                n2sq = sum(y.^2,1);
                n2 = size(y,2);
            end;
            K = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq - 2*x'*y;
            K(K<0) = 0;
            K = exp(-K/this.Gamma);
        end
                
        function dx = evaluateD1(this, x)
            % Method for first derivative evaluation
            dx = -2*x/this.Gamma .* exp(-x.^2/this.Gamma);
        end
        
        function ddx = evaluateD2(this, x)
            % Method for second derivative evaluation
            ddx = (2/this.Gamma) * (2*x.^2/this.Gamma-1) .* exp(-x.^2/this.Gamma);
        end
        
        function Kx = evaluateScalar(this, x)
            % Implements the required method from the IRotationInvariant
            % interface
            Kx = exp(-x.^2/this.Gamma);
        end
        
        function set.Gamma(this, value)
            if ~isposrealscalar(value)
                error('Only positive scalar values allowed for Gamma.');
            end
            this.Gamma = value;
            % Adjust the BellFunctions' x0 value
            this.x0 = sqrt(value/2); %#ok
            this.PenaltyFactor = 1/value; %#ok
        end
        
        function g = setGammaForDistance(this, dist, ep)
            % Computes the `\gamma` value for which the Gaussian is smaller
            % than `\epsilon` in a distance of dist. Returns the computed
            % value AND sets the kernel's Gamma property to this value.
            %
            % Parameters:
            % dist: The target distance at which the gaussian is smaller
            % than ep
            % ep: The `\epsilon` value. If not given, `\epsilon`=eps
            % (machine precision) is assumed.
            %
            % Return values:
            % g: The computed gamma
            if nargin == 2
                ep = eps;
            end
            g = -(dist^2)/log(ep);
            this.Gamma = g;
        end
    end
    
end

