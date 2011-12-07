classdef AMassMatrix < KerMorObject & dscomponents.IProjectable
% AMassMatrix: 
%
%
%
% @author Daniel Wirtz @date 2011-12-06
%
% @new{0,6,dw,2011-12-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
%
%
% @todo return P,Q if M is sparse and include a flag for sparsity of M
% (advantage? is done in matlabs ODE23 solvers)
    
    methods(Abstract)
        M = evaluate(this, t);
        
        [L,U,Q,P] = getLU(this, t);
    end
    
end