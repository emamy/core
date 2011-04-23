classdef KernelLS < KerMorObject & approx.IKernelCoeffComp
    %KERNELLS Least-Squares kernel regression ("Rigde Regression")
    %   Since the systems can be considerably large, the pcg solver is used
    %   instead of plain inversion.
    %
    % @author Daniel Wirtz 21.03.2010
    
    properties
        % The kernel matrix to use für LS regression
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
                % @TODO: why does pcg not work here? matrix is symmetric!
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
        
    end
    
end

