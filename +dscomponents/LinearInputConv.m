classdef LinearInputConv < dscomponents.AInputConv
    % Simple linear (=matrix) input conversion
    
    properties(Access=private)
        % The target matrix
        B;
    end
    
    methods
        function this = LinearInputConv(B)
            this.B = B;
        end
        
        function res = evaluate(this, t, mu)%#ok
            res = this.B;
        end
        
        function proj = project(this, V, W)%#ok
            proj = dscomponents.LinearInputConv(W' * this.B);
        end
        
        function copy = clone(this)
            copy = dscomponents.LinearInputConv(this.B);
        end
    end
    
end

