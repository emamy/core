classdef KernelLS < KerMorObject & approx.IKernelCoeffComp
    %KERNELLS Least-Squares kernel regression ("Rigde Regression")
    %   Since the systems can be considerably large, the pcg solver is used
    %   instead of plain inversion.
    %
    % @author Daniel Wirtz 21.03.2010
    %
    % @change{0,4,dw,2011-05-03} Removed `b` offset terms as no longer used in kernel expansions.
    %
    % @todo add unit test
    
    % @change{0,4,sa,2011-05-07} Implemented Setter for the properties G and Epsilon
    
    properties(SetObservable)
        % The kernel matrix to use für LS regression
        %
        % @propclass{data} Needed for LS to work in the first place.
        K;
        
        % The regularization term weight
        %
        % @propclass{important} If the data is badly conditioned increasing lambda can help.
        %
        % @default 1
        lambda = 1;
        
        % Maximum iteration count for the conjugate gradient method
        %
        % @propclass{alglimit} Limits the iterations for the pcg method.
        %
        % See also: pcg, bicg
        CGMaxIt = [];
        
        % Tolerance for cg method
        %
        % @propclass{critical} The error tolerance for the pcg method.
        %
        % See also: pcg, bicg
        CGTol = 1e-6;
        
        % Maximum dimension of function for which direct inversion is
        % performed instead of pcg
        %
        % @propclass{optional} Direct inversion is more precise for small sizes. 
        %
        % @default 2000
        MaxStraightInvDim = 2000;
    end
    
    methods
        
        function this = KernelLS
            this = this@KerMorObject;
            this.registerProps('K','lambda','CGMaxIt','CGTol','MaxStraightInvDim');
        end
        
        function a = regress(this,fx)
            
            % Ensure fxi is a column vector
            fx = reshape(fx,size(this.K,1),[]);
            
            y = this.K'*fx;
            M = this.K'*this.K + this.lambda*eye(size(this.K));
            
            if length(y) <= this.MaxStraightInvDim
                a = M\y;
            else
                % @TODO: why doepcgs pcg not work here? matrix is symmetric!
                [a, flag] = bicg(M,y,this.CGTol, this.CGMaxIt);
            end
        end
        
        %% approx.IKernelCoeffComp interface members
        function init(this, K)
            this.K = K;
        end
        
        function [ai, svidx] = computeKernelCoefficients(this, yi)
            ai = this.regress(yi);
            svidx = [];
        end
        
        %% setters
        
        function set.K(this, value)
            if ~isa(value, 'double')
                error('value must be a double matrix');
            end
            % Make matrix symmetric (can be false due to rounding errors)
            this.K = .5*(value + value');
        end
        
        function set.lambda(this, value)
            if ~isposintscalar(value)
                error('value must be a positive integer scalar');
            end
            this.lambda = value;
        end
        
        function set.CGMaxIt(this, value)
            if ~isposintscalar(value)
                error('value must be a positive integer scalar');                
            end
            this.CGMaxIt = value;
        end
        
        function set.CGTol(this, value)
            if ~isposrealscalar(value)
                error('value must be a positive real scalar.');                
            end
            this.CGTol = value;
        end
        
        function set.MaxStraightInvDim(this, value)
            if ~isposintscalar(value)
                error('value must be a positive integer scalar');                
            end
            this.MaxStraightInvDim = value;
        end
    end
    
end
