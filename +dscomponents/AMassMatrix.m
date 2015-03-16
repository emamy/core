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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

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
    
    properties(SetAccess=private)
        AlgebraicEquationDims = [];
    end
    
    methods
        function this = AMassMatrix(algdims, varargin)
            this = this@KerMorObject;
            if nargin == 1
                this.AlgebraicEquationDims = algdims;
            end
        end
        
        function copy = clone(this, copy)
            copy = clone@general.AProjectable(this, copy);
            copy.TimeDependent = this.TimeDependent;
            copy.SparsityPattern = this.SparsityPattern;
            copy.AlgebraicEquationDims = this.AlgebraicEquationDims;
        end
    end
    
    methods(Abstract)
        M = evaluate(this, t, mu);
        
        [L,U,Q,P] = getLU(this, t, mu);
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if nargin == 2
                obj.TimeDependent = from.TimeDependent;
                obj.SparsityPattern = from.SparsityPattern;
                obj = loadobj@KerMorObject(obj, from);
                obj = loadobj@general.AProjectable(obj, from);
            else
                obj = loadobj@KerMorObject(obj);
                obj = loadobj@general.AProjectable(obj);
            end
        end
    end
    
end