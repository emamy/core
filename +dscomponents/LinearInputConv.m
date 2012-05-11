classdef LinearInputConv < dscomponents.AInputConv
    % Simple linear (=matrix) input conversion
    
    properties(SetAccess=private)
        % The target matrix
        B;
    end
    
    methods
        function this = LinearInputConv(B)
            this.B = B;
        end
        
        function res = evaluate(this, ~, ~)
            res = this.B;
        end
        
        function proj = project(this, ~, W)
            proj = dscomponents.LinearInputConv(W' * this.B);
        end
        
        function copy = clone(this)
            copy = dscomponents.LinearInputConv(this.B);
        end
        
        function prod = mtimes(this, other)
            prod = dscomponents.LinearInputConv(this.B * other.B);
        end
        
        function diff = minus(this, other)
            diff = dscomponents.LinearInputConv(this.B - other.B);
        end
        
        function diff = plus(this, other)
            diff = dscomponents.LinearInputConv(this.B + other.B);
        end
        
        function transp = ctranspose(this)
            transp = dscomponents.LinearInputConv(this.B');
        end
    end
    
end

