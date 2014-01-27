classdef Constant < error.alpha.Base
% Constant: Constant alpha terms
%
% @author Daniel Wirtz @date 2011-07-04
%
% @new{0,5,dw,2011-07-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Access=private)
        M2 = [];
        M3 = [];
    end
    
    methods
        
        function this = Constant(model)
            this = this@error.alpha.Base(model);
        end
        
        function inputOfflineComputations(this, rm, M)
            % Performs the offline stage for the error estimators regarding
            % the inputs.
            %
            % Parameters:
            % rm: The reduced model @type models.ReducedModel
            % M: The projected coefficient matrix `\vM_{\alpha} -
            % \vV\vW^T\vM_{\alpha})` @type matrix<double>
            fm = rm.FullModel;
            if ~isempty(fm.System.B)
                try
                    B = fm.System.B.evaluate([],[]);
                catch ME%#ok
                    B = fm.System.B.evaluate(0,rm.getRandomParam);
                    warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                end
            
                if ~isempty(rm.V)
                    B2 = B - rm.V*(rm.W'*B);
                    this.M2 = M'*(rm.G*B2);
                    this.M3 = B2'*(rm.G*B2);
                    clear B2;
                else
                    n = size(this.M1,2);
                    b = size(B,2);
                    this.M2 = zeros(n,b);
                    this.M3 = zeros(b,b);
                end
            end
        end
        
        function a = getAlpha(this, phi, ut, t, mu)%#ok
            % Computes the alpha term for the error estimator
            %
            % Parameters:
            % phi: The kernel vector `\K(x,x_i)` @type colvec
            % ut: The evaluation of the current input `u(t)` at time `t`
            % @type double
            % t: The current time `t` @type double
            % mu: The current parameter `\mu` @type colvec
            a = phi*this.M1*phi';
            if ~isempty(ut) % An input function u is set
                a = a + phi*this.M2*ut + ut'*this.M3*ut;
            end
            a = sqrt(abs(a));
            %a = sqrt(max(a,0));
        end
    end
    
end