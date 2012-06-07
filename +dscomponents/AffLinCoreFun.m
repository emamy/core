classdef AffLinCoreFun < dscomponents.ACoreFun & general.AffParamMatrix ...
        & dscomponents.IGlobalLipschitz
%Simple affine-linear core function "f" for a dynamical system.
% 
% Simply wraps an affine-linear function into the ACoreFun interface to
% enable use of simple affine-linear functions as core function. At
% projection, each summand matrix is base changed into the basis given
% by V.
%
% @author Daniel Wirtz @date 15.03.2010
%
% @change{0,6,dw,2012-02-02} Removed the former offset term `b` as it can be modeled via the
% input source `B(t,\mu)u(t)` with constant `u \equiv 1`.
%
% @new{0,6,dw,2011} Added an optional offset term `b` to the AffLinCoreFun
% to enable affine-linear affine-parametric core functions.
%
% @change{0,5,dw,2011-07-07} Updated this class to use the new general.AffParamMatrix class inside.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetObservable)
        % Export setting. Java class name for JKerMor model export
        %
        % Set this value to the class inside your JKerMor source that
        % implements the IAffineCoefficients interface for this core
        % function.
        %
        % In order for this to work the coefficient functions must equal
        % the AffParamMatrix' coefficient functions both mathematically and
        % in the order of entry.
        %
        % @propclass{data} Set only if the model is intended for
        % JKerMor export.
        %
        % @type char @default ''
        CoeffClass = '';
    end
   
    methods
        function this = AffLinCoreFun
            % Creates a new instance of the AffLinCoreFun.
            this = this@general.AffParamMatrix;
            this = this@dscomponents.ACoreFun;
            this.CustomProjection = true;
        end
        
        function fx = evaluate(this, x, t, mu)
            fx = this.compose(t, mu)*x;
        end
        
        function fx = evaluateCoreFun(varargin)%#ok
            % This method will never be called as evaluate is overridden
            % directly for performance.
            % this is possible as AffLinCoreFuns have both CustomProjection
            % and MultiArgumentEvaluations.
            error('This method should never be called.');
        end
        
        function J = getStateJacobian(this, ~, t, mu)
            % Overrides the default jacobian finite difference
            % implementation.
            J = this.compose(t, mu);
        end
        
        function c = getGlobalLipschitz(this, t, mu)
            % Implementation of the interface method from IGlobalLipschitz.
            error('need to update.');
            a = this.AffParamMatrix;
            c = 0;
            for idx=1:length(this.Coefficients)
                cfun = this.Coefficients{idx};
                c = c + abs(cfun(t,mu)) * norm(a.Matrices(:,:,idx));
            end
        end
        
        function proj = project(this, V, W)
            proj = this.clone;
            proj = project@dscomponents.ACoreFun(this, V, W, proj);
            proj = project@general.AffParamMatrix(this, V, W, proj);
            proj.JSparsityPattern = [];
        end
        
        function addMatrix(this, coeff_fcn, mat)
            % Update the XDim as first matrix is added
            if isempty(this.XDim)
                this.XDim = size(mat,2);
            end
            
            addMatrix@general.AffParamMatrix(this, coeff_fcn, mat);
            
            % Compute sparsity pattern
            if issparse(mat)
                if ~isempty(this.JSparsityPattern)
                    this.JSparsityPattern = mat ~= 0 && this.JSparsityPattern;
                else
                    this.JSparsityPattern = mat ~= 0;
                end
            else
                this.JSparsityPattern = [];
            end
            
            mu = ones(1,100);
            this.TimeDependent = ~all(this.cfun(0,mu) == this.cfun(Inf,mu));
            fprintf('AffLinCoreFun: Guessed time-dependency to %d.\n',this.TimeDependent);
        end
        
        function prod = mtimes(this, other)
            prod = mtimes@general.AffParamMatrix(this, other);
            prod.postprocess;
        end
        
        function diff = minus(this, other)
            diff = minus@general.AffParamMatrix(this, other);
            diff.postprocess;
        end
        
        function sum = plus(this, other)
            sum = plus@general.AffParamMatrix(this, other);
            sum.postprocess;
        end
        
        function transp = ctranspose(this)
            transp = ctranspose@general.AffParamMatrix(this);
            transp.postprocess;
        end
    end
    
    methods(Access=private)
        function postprocess(this)
            if this.N > 0
                this.XDim = this.dims(2);
                this.JSparsityPattern = [];
            end
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = dscomponents.AffLinCoreFun;
            copy = clone@general.AffParamMatrix(this, copy);
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.CoeffClass = this.CoeffClass;
        end
    end
end

