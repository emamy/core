classdef ConstMassMatrix < dscomponents.AMassMatrix
% ConstMassMatrix: 
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
    
    properties(SetAccess=private)
        M;
    end
    
    properties(Access=private)
        l;
        u;
        q;
        p;
    end
    
    methods
        
        function this = ConstMassMatrix(M)
            this.M = M;
            this.TimeDependent = false;
            [this.l, this.u, this.q, this.p] = lu(M);
            % Compute sparsity pattern straight away
            [i, j] = find(M);
            this.SparsityPattern = sparse(i,j,ones(size(i)));
        end
        
        function M = evaluate(this, ~)
            M = this.M;
        end
        
        function [L, U, q, p] = getLU(this, ~, ~)
            L = this.l;
            U = this.u;
            q = this.q;
            p = this.p;
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected.M = W'*(this.M*V);
        end
        
        function copy = clone(this)
            copy = dscomponents.ConstMassMatrix(this.M);
        end
    end
    
end