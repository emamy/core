classdef AffLinInputConv < general.AffParamMatrix & dscomponents.AInputConv
% AffLinInputConv: Affine parametric input conversion matrix `B(t,\mu)`
%
% This matrix has the structure `B(t,\mu) = \sum\limits_{i=0}^Q \theta_i(t,\mu)B_i`.
%
% Basically extends the general.AffParamMatrix and wraps it into the dscomponents.AInputConv
% interface.
%
% See also: general.AffParamMatrix
%
% @new{0,5,dw,2011-07-04} Added this class.
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
    
    properties(Access=private)
        cachedB;
    end
    
    methods
        function B = evaluate(this, t, mu)
            % Evaluates the input conversion matrix.
            %
            % For this case, it is simply calling compose of the superclass AffParamMatrix.
            %
            % Parameters:
            % t: The current time `t`
            % mu: The current parameter vector `\mu`
            %
            % Return values:
            % B: The affine parametric matrix `B(t,\mu)`.
            if this.TimeDependent
                B = this.compose(t, mu);
            else
                B = this.cachedB;
            end
        end
        
        function prepareSimulation(this, mu)
            this.cachedB = [];
            if ~isempty(mu) && ~this.TimeDependent
                this.cachedB = this.compose(0, mu);
            end
        end
        
        function projected = project(this, V, W)%#ok
            % Projects the affine parametric input conversion matrix B into the subspace spanned by
            % `V,W`.
            
            % Uses the overridden operators in AffParamMatrix to create a copy.
            projected = W'*this;
            
            % Dont need to store V,W for disk space reasons
            %projected = project@general.AProjectable(this, V, W, projected);
        end
        
        function copy = clone(this)
            copy = dscomponents.AffLinInputConv;
            copy = clone@general.AffParamMatrix(this, copy);
            copy.CoeffClass = this.CoeffClass;
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, varargin)
            obj = loadobj@general.AffParamMatrix(obj, varargin{:});
            obj = loadobj@dscomponents.AInputConv(obj, varargin{:});
        end
    end
    
end

