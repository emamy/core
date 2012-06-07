classdef AMassMatrix < KerMorObject & general.AProjectable
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

    properties(SetAccess=protected)
        % Flag that indicates time-dependency of the Mass Matrix.
        %
        % Set in subclasses to mirror correct behaviour.
        %
        % @type logical @default false
        TimeDependent = false;
        
        % The sparsity pattern for the mass matrix
        %
        % Set this value in subclasses if the mass matrix has a specific
        % sparsity pattern. This can speed up computations for large
        % systems with mass matrix considerably.
        %
        % @type sparsematrix @default []
        SparsityPattern = [];
    end
    
    methods
        function copy = clone(this, copy)
            copy = clone@general.AProjectable(this, copy);
            copy.TimeDependent = this.TimeDependent;
            copy.SparsityPattern = this.SparsityPattern;
        end
    end
    
    methods(Abstract)
        M = evaluate(this, t, mu);
        
        [L,U,Q,P] = getLU(this, t, mu);
    end
    
end