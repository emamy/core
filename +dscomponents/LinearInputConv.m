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
        
        function proj = project(this, V, W)
            proj = dscomponents.LinearInputConv(W' * this.B);
            proj = project@general.AProjectable(this, V, W, proj);
        end
        
        function copy = clone(this)
            copy = dscomponents.LinearInputConv(this.B);
        end
        
        function prod = mtimes(this, other)
            if isa(this,'dscomponents.LinearInputConv')
                if isa(other,'dscomponents.LinearInputConv')
                    prod = dscomponents.LinearInputConv(this.B * other.B);
                else
                    prod = dscomponents.LinearInputConv(this.B * other);
                end
            elseif isa(other,'dscomponents.LinearInputConv')
                prod = dscomponents.LinearInputConv(this * other.B);
            end
        end
        
        function diff = minus(this, other)
            if isa(this,'dscomponents.LinearInputConv')
                if isa(other,'dscomponents.LinearInputConv')
                    diff = dscomponents.LinearInputConv(this.B - other.B);
                else
                    diff = dscomponents.LinearInputConv(this.B - other);
                end
            elseif isa(other,'dscomponents.LinearInputConv')
                diff = dscomponents.LinearInputConv(this - other.B);
            end
        end
        
        function diff = plus(this, other)
            if isa(this,'dscomponents.LinearInputConv')
                if isa(other,'dscomponents.LinearInputConv')
                    diff = dscomponents.LinearInputConv(this.B + other.B);
                else
                    diff = dscomponents.LinearInputConv(this.B + other);
                end
            elseif isa(other,'dscomponents.LinearInputConv')
                diff = dscomponents.LinearInputConv(this + other.B);
            end
        end
        
        function transp = ctranspose(this)
            transp = dscomponents.LinearInputConv(this.B');
        end
    end
    
end

