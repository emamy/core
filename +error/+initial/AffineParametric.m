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
        m;
    end
    
    methods
        function prepareForReducedModel(this, rm)
            % Operates on affine-parametric matrices, see general.AffParamMatrix
            this.m = rm;
            % Consider scaling, the state space variables are scaled if set, so need the
            % correct scaled x0 value here.
%             s = model.System.StateScaling;
%             if s ~= 1
%                 warning('KerMor:unchecked','Functionality not yet checked for scaled systems.');
%                 if isscalar(s)
%                     dim = size(model.System.x0.getMatrix(1),2);
%                     s(1:dim,1) = s;
%                 else
%                     dim = length(s);
%                 end
%                 x0 = model.System.x0 * spdiags(1./s,0,dim,dim);
%             else
                x0 = rm.FullModel.System.x0;
%             end
            this.af = x0 - rm.V*(rm.W'*x0);
        end
        
        function e0 = getE0(this, mu)
            e0 = this.af.compose(mu);
            e0 = Norm.LG(e0,this.m.G);
        end
    end
    
end