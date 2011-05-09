classdef KernelLS < KerMorObject & approx.IKernelCoeffComp
    %KERNELLS Least-Squares kernel regression ("Rigde Regression")
    %   Since the systems can be considerably large, the pcg solver is used
    %   instead of plain inversion.
    %
    % @author Daniel Wirtz 21.03.2010
    %
    % @change{0,3,sa,2011-05-07} Implemented Setter for the properties G
    % and Epsilon
    
    properties
        % The kernel matrix to use f�r LS regression
        K;
        
        % The regularization term weight
        lambda=1;
        
        % Maximum iteration count for the conjugate gradient method
        % 
        % See also: pcg, bicg
        CGMaxIt = [];
        
        % Tolerance for cg method
        % 
        % See also: pcg, bicg
        CGTol = 1e-6;
        
        % Maximum dimension of function for which direct inversion is
        % performed instead of pcg
        %
        % Defaults to 2000.
        MaxStraightInvDim = 2000;
    end
    
    methods
        
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
        
        function [ai, b, svidx] = computeKernelCoefficients(this, yi)
            ai = this.regress(yi);
            b = 0;
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
