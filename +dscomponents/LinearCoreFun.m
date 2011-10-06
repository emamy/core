classdef LinearCoreFun < dscomponents.ACoreFun
    % Linear core function for state space models (no time or parameter)
    
    properties(SetAccess=protected)
        % The systems core function matrix
        A;
    end
    
    methods
        function this = LinearCoreFun(A)
            if nargin == 1
                this.A = A;
            end
            this.CustomProjection = true;
            this.MultiArgumentEvaluations = true;
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
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            fx = this.A*x;
        end
    end
    
end

