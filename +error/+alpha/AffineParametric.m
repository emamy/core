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
        
        function this = AffineParametric(rmodel)
            this = this@error.alpha.Base(rmodel);
        end
        
        function inputOfflineComputations(this, rmodel, M)
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                fm = rmodel.FullModel;
                % Please see the AffParamMatrix overridden operators to understand what's going on
                % here :-)
                B = fm.System.B - rmodel.V*(rmodel.W'*fm.System.B);
                this.mb = (M'*rmodel.GScaled)*B;
                this.bb = B'*(rmodel.GScaled*B);
            end
        end
        
        function a = getAlpha(this, phi, ut, t, mu)
            a = phi*this.M1*phi';
            if ~isempty(ut) % An input function u is set
                a = a + phi*this.mb.compose(t,mu)*ut + ut'*this.bb.compose(t,mu)*ut;
            end
            a = sqrt(abs(a));
        end
       
    end
    
end