classdef AffineParametric < error.alpha.Base
    % AffineParametric:
    %
    %
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
        mb = [];
        bb = [];
    end
    
    methods
        
        function this = AffineParametric(rm)
            this = this@error.alpha.Base(rm);
        end
        
        function inputOfflineComputations(this, rm, M)
            % Performs the offline stage for the error estimators regarding
            % the inputs.
            %
            % Parameters:
            % rm: The reduced model @type models.ReducedModel
            % M: The projected coefficient matrix `\vM_{\alpha} - \vV\vW^T\vM_{\alpha})` @type
            % matrix<double>
            fm = rm.FullModel;
            if ~isempty(rm.V)
                % For details see the AffParamMatrix overridden operators
                B = fm.System.B - rm.V*(rm.W'*fm.System.B);
                this.mb = (M'*rm.G)*B;
                this.bb = B'*(rm.G*B);
            end
        end
        
        function a = getAlpha(this, phi, ut, t, mu)
            % Computes the `\alpha(x, t,\mu)` value of the error
            % estimator.
            %
            % Parameters:
            % phi: The kernel vector `\K(x,x_i)` @type colvec
            % ut: The evaluation of the current input `\vu(t)` at time `t` @type colvec<double>
            % t: The current time `t` @type double
            % mu: The current parameter `\vmu` @type colvec<double>
            a = phi*this.M1*phi';
            if ~isempty(ut) % An input function u is set
                a = a + phi*this.mb.compose(t,mu)*ut + ut'*this.bb.compose(t,mu)*ut;
            end
            a = sqrt(abs(a));
        end
       
    end
    
end