classdef LinearCoreFun < dscomponents.ACoreFun
    % Linear core function for state space models (no time or parameter)
    %
    % @new{0,6,dw,2011-12-7} Implemented the IJacobian interface.
    
    properties(SetAccess=protected)
        % The systems core function matrix
        A;
    end
    
    methods
        function this = LinearCoreFun(A)
            if nargin > 0
                this.A = A;
                % If not set here, subclasses must take care to set xDim
                % whenever A changes.
                this.fDim = size(A,1);
                this.xDim = size(A,2);
                % Compute sparsity pattern
                if issparse(A)
                    this.JSparsityPattern = A ~= 0;
                end
            end
            this.CustomProjection = true;
            this.MultiArgumentEvaluations = true;
            this.TimeDependent = false;
        end
        
        function projected = project(this, V, W, projected)
            if nargin < 4
                projected = this.clone;
            end
            projected.A = W'*(this.A*V);
        end
        
        function copy = clone(this)
            copy = clone@dscomponents.ACoreFun(this, dscomponents.LinearCoreFun);
            copy.A = this.A;
        end
        
        function fx = evaluate(this, x, ~, ~)
            fx = this.A*x;
        end
        
        function fx = evaluateCoreFun(varargin)%#ok
            % This method will never be called as evaluate is overridden
            % directly for performance.
            % this is possible as LinearCoreFuns have both CustomProjection
            % and MultiArgumentEvaluations.
            error('This method should never be called.');
        end
        
        function J = getStateJacobian(this, ~, ~, ~)
            J = this.A;
        end
        
        function prod = mtimes(this, other)
            if isa(this,'dscomponents.LinearCoreFun')
                if isa(other,'dscomponents.LinearCoreFun')
                    prod = dscomponents.LinearCoreFun(this.A * other.A);
                else
                    prod = dscomponents.LinearCoreFun(this.A * other);
                end
            elseif isa(other,'dscomponents.LinearCoreFun')
                prod = dscomponents.LinearCoreFun(this * other.A);
            end
        end
        
        function diff = minus(this, other)
            diff = dscomponents.LinearCoreFun(this.A - other.A);
        end
        
        function diff = plus(this, other)
            diff = dscomponents.LinearCoreFun(this.A + other.A);
        end
        
        function transp = ctranspose(this)
            transp = dscomponents.LinearCoreFun(this.A');
        end
    end
end

