classdef InvMultiquadrics < kernels.ARBFKernel
    %InvMULTIQUADRICS Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @todo check for what params the inv multiquadrics are a bell function
    % and implement bell function methods!
    
    properties(SetObservable)
        % The negative exponent
        %
        % @propclass{critical} Greatly influences the kernels behaviour.
        beta=-1;
        
        % The offset value
        %
        % @propclass{critical} Greatly influences the kernels behaviour.
        c=1;
    end
    
    methods
        
        function this =  InvMultiquadrics(beta, c)
            % Constructor offering the possibility to initialize kernel
            % specifics at creation time.
            
            this = this@kernels.ARBFKernel;
            
            % Register before processing arguments, because if set that's a custom user option.
            this.registerProps('beta','c');
            
            if nargin > 0
                this.beta = beta;
                if nargin > 1
                    this.c = c;
                end
            end
        end
        
        function c = getGlobalLipschitz(this)%#ok
            % @todo implement
            error('Not implemented yet');
        end
        
        function Nabla = getNabla(this, x, y)%#ok
            % @todo implement
            error('Not implemented yet');
        end
        
        function K = evaluate(this, x, y)
            % Evaluates the inverse multiquadrics kernel.
            %
            % If `y_j` is set, the dimensions of `x_i` and `y_j` must be equal for all `i,j`.
            %
            % Parameters:
            % x: First set `x_i \in \R^d` of `n` vectors @type matrix<double>
            % y: Second set `y_j \in \R^d` of `m` vectors. If y is empty `y_i = x_i` and `n=m`
            % is assumed. @type matrix<double>
            %
            % Return values:
            % K: An evaluation matrix `K \in \R^{n\times m}` of the evaluated multiquadrics
            % with entries `K_{i,j} = 1/(c^2+\epsilon\norm{x_i-y_j}{G}^2)`.
            K = (this.c^2 + this.epsilon^2*this.getSqDiffNorm(x,y)).^this.beta;
        end
        
        function Ks = evaluateScalar(this, r)
            Ks = (this.c^2 + (this.epsilon*r).^2).^this.beta;
        end
    end
    
end

