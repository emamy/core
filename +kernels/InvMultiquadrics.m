classdef InvMultiquadrics < kernels.BaseKernel & kernels.ARotationInvariant
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
            
            this = this@kernels.BaseKernel;
            this = this@kernels.ARotationInvariant;
            
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
        
%         function K = evaluate(this, x, y)
%             % Original source from Bernard Haasdonk / KerMet
%             n1sq = sum(x.^2,1);
%             n1 = size(x,2);
%             
%             if nargin == 2;
%                 n2sq = n1sq;
%                 n2 = n1;
%                 y = x;
%             else
%                 n2sq = sum(y.^2,1);
%                 n2 = size(y,2);
%             end;
%             r = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq - 2*x'*y;
%             K = (this.c^2 + r.^2).^this.beta;
%         end
        
        function Ks = evaluateScalar(this, x)
            Ks = (this.c^2 + x.^2).^this.beta;
        end
    end
    
end

