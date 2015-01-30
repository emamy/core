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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
    
    properties(Transient, Access=private)
        cachedA;
    end
   
    methods
        function this = AffLinCoreFun(sys)
            % Creates a new instance of the AffLinCoreFun.
            this = this@general.AffParamMatrix;
            this = this@dscomponents.ACoreFun(sys);
            this.CustomProjection = true;
            this.TimeDependent = false;
            this.cachedA = [];
        end
        
        function fx = evaluate(this, x, t)
            % Evaluates affine-linear core function by matrix-vector
            % multiplication.
            if this.TimeDependent
                fx = this.compose(t, this.mu)*x;
            else
                fx = this.cachedA*x;
            end
        end
        
        function fx = evaluateMulti(this, x, t, mu)
            if (~this.TimeDependent || numel(t) == 1) && size(mu,2) == 1
                fx = this.compose(t, mu)*x;
            else
                % Multi-Argument case
                fx = zeros(size(x));
                for k = 1:size(x,2)
                    fx(:,k) = this.compose(t(k), mu(:,k))*x(:,k);
                end
            end
        end
        
        function prepareSimulation(this, mu) 
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            this.cachedA = [];
            if ~isempty(mu) && ~this.TimeDependent
                this.cachedA = this.compose(0, mu);
            end
        end
        
        function fx = evaluateCoreFun(varargin)%#ok
            % This method will never be called as evaluate is overridden
            % directly for performance.
            % this is possible as AffLinCoreFuns have both CustomProjection
            % and MultiArgumentEvaluations.
            error('This method should never be called.');
        end
        
        function J = getStateJacobian(this, ~, t)
            % Overrides the default jacobian finite difference
            % implementation. Jacobian of linear operator is the operator
            % itself, i.e. in this case the linear combination of the matrices.
            if this.TimeDependent
                J = this.compose(t, this.mu);
            else
                J = this.cachedA;
            end
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
            % Adds a new matrix to the affine-linear core function.
            %
            % Parameters:
            % coeff_fcn: A string describing the `i`-th coefficient function's evaluation as if
            % entered into a function handle interior or a function. By convention, the passed arguments are named
            % 't' and 'mu'. @type string
            % mat: The corresponding matrix `A_i` @type matrix<double>
            %
            % See also: general.AffParramMatrix
            
            % Update the xDim as first matrix is added
            if isempty(this.xDim)
                this.xDim = size(mat,2);
            end
            if isempty(this.fDim)
                this.fDim = size(mat,1);
            end
            
            addMatrix@general.AffParamMatrix(this, coeff_fcn, mat);
            
            % Compute sparsity pattern
            if issparse(mat)
                if ~isempty(this.JSparsityPattern)
                    this.JSparsityPattern = mat ~= 0 | this.JSparsityPattern;
                else
                    this.JSparsityPattern = mat ~= 0;
                end
            else
                this.JSparsityPattern = [];
            end
            
            mu = ones(1,100);
            if ~this.TimeDependent
                this.TimeDependent = ~all(this.cfun(0,mu) == this.cfun(Inf,mu));
                if KerMor.App.Verbose > 1
                    fprintf('AffLinCoreFun: Guessed time-dependency to %d.\n',this.TimeDependent);
                end
            end
        end
        
        function prod = mtimes(this, other)
            % Implements the default multiplication method.
            prod = mtimes@general.AffParamMatrix(this, other);
            prod.postprocess;
        end
        
        function diff = minus(this, other)
            % Implements the default subtraction method.
            diff = minus@general.AffParamMatrix(this, other);
            diff.postprocess;
        end
        
        function sum = plus(this, other)
            % Implements the default addition method.
            sum = plus@general.AffParamMatrix(this, other);
            sum.postprocess;
        end
        
        function transp = ctranspose(this)
            % Implements the default transposition method.
            transp = ctranspose@general.AffParamMatrix(this);
            transp.postprocess;
        end
    end
    
    methods(Access=private)
        function postprocess(this)
            % If there is a affine-linear core function, update x-dimension
            % and reset Jacobian SparsityPattern after operation on
            % affine-linear core function was performed.
            if this.N > 0
                this.xDim = this.dims(2);
                this.JSparsityPattern = [];
            end
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = dscomponents.AffLinCoreFun(this.System);
            copy = clone@general.AffParamMatrix(this, copy);
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.CoeffClass = this.CoeffClass;
        end
    end
end

