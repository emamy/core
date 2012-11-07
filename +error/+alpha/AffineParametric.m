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
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Access=private)
        mb = [];
        bb = [];
    end
    
    methods
        
        function this = AffineParametric(fm)
            this = this@error.alpha.Base(fm);
        end
        
        function inputOfflineComputations(this, fm, M)
            % Performs the offline stage for the error estimators regarding
            % the inputs.
            %
            % Parameters:
            % fm: The full model @type models.BaseFullModel
            % M: The projected coefficient matrix `\vM_{\alpha} - \vV\vW^T\vM_{\alpha})` @type
            % matrix<double>
            if ~isempty(fm.Data.V) && ~isempty(fm.Data.W)
                % Please see the AffParamMatrix overridden operators to understand what's going on
                % here :-)
                B = fm.System.B - fm.Data.V*(fm.Data.W'*fm.System.B);
                this.mb = (M'*fm.GScaled)*B;
                this.bb = B'*(fm.GScaled*B);
            end
        end
        
        function a = getAlpha(this, phi, ut, t, mu)
            % Computes the `\alpha(x, t,\mu)` value of the error
            % estimator.
            %
            % Parameters:
            % phi: The kernel vector `\Phi(x,x_i)` @type colvec
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