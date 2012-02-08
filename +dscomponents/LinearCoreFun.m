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
            end
            this.CustomProjection = true;
            this.MultiArgumentEvaluations = true;
            this.TimeDependent = false;
            this.CustomJacobian = true;
        end
        
        function projected = project(this, V, W, projected)
            if nargin < 4
                projected = this.clone;
            end
            projected.A = W'*(this.A*V);
        end
        
        function copy = clone(this, copy)
            if nargin == 1
                copy = dscomponents.LinearCoreFun;
            end
            copy.A = this.A;
            copy = clone@dscomponents.ACoreFun(this, copy);
        end
        
        function fx = evaluateCoreFun(this, x, ~, ~)
            fx = this.A*x;
        end
        
        function J = getStateJacobian(this, ~, ~, ~)
            J = this.A;
        end
    end
    
end

