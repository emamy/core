classdef LinearCoreFun < dscomponents.ACoreFun & dscomponents.IJacobian
    % Linear core function for state space models (no time or parameter)
    %
    % @new{0,6,dw,2011-12-7} Implemented the IJacobian interface.
    
    properties(SetAccess=protected)
        % The systems core function matrix
        A;
        
        % The linear core function's offset vector
        b = [];
    end
    
    methods
        function this = LinearCoreFun(A, b)
            if nargin > 0
                this.A = A;
                if nargin > 1
                    this.b = b;
                end
            end
            this.CustomProjection = true;
            this.MultiArgumentEvaluations = true;
        end
        
        function projected = project(this, V, W, projected)
            if nargin < 4
                projected = this.clone;
            end
            projected.A = W'*(this.A*V);
            if ~isempty(this.b)
                projected.b = W'*this.b;
            end
        end
        
        function copy = clone(this, copy)
            if nargin == 1
                copy = dscomponents.LinearCoreFun;
            end
            copy.A = this.A;
            copy.b = this.b;
            copy = clone@dscomponents.ACoreFun(this, copy);
        end
        
        function fx = evaluateCoreFun(this, x, ~, ~)
            fx = this.A*x;
            if ~isempty(this.b)
                fx = fx + this.b;
            end
        end
        
        function J = getStateJacobian(this, ~, ~, ~)
            J = this.A;
        end
    end
    
end

