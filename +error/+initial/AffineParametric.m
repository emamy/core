classdef AffineParametric < error.initial.Base
% AffineParametric: 
%
% Computes `\tilde{x_i} = (I_d - VW^t)x_i` beforehand and 
% `\no{\sum\limits_{i=0}^{Q_0}\tilde{x_i}}` at online stage.
% Due to that it is not really completely offline-online separated, but as the coefficient functions
% may be negative this is the better estimation.
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
        af;
        rm;
    end
    
    methods
        function this = AffineParametric(rmodel)
            this.af = general.AffParamMatrix;
            this.rm = rmodel;
            x0 = rmodel.FullModel.System.x0;
            % Operate on affine-parametric matrices, see general.AffParamMatrix
            this.af = x0 - rmodel.V*(rmodel.W'*x0);
        end
        
        function e0 = getE0(this, mu)
            e0 = this.af.compose(mu);
            e0 = sqrt(e0'*this.rm.GScaled*e0);
        end
    end
    
end