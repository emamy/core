classdef DoxyGenReferenceClass < testing.MUnit & handle
    %DOXYGENREFERENCECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private, GetAccess=protected)
       
       % Some comment!
       % @type testing.MUnit
       SomeClass = testint.MUnit;
    end
    
    properties
       % Some comment
       % @type int
       SomeProp = 0;
    end
    
    properties(Dependent)
       % equals someprop times five. 
       SomeDepProp;
    end
    
    methods(Sealed,Access=private)
        
    end
    
    methods
        % Some function 
        function set.SomeProp(this, value)
            if value > 0
                this.SomeProp = value;
            else
                error('buh');
            end
        end
        
        % Some redundant method.
        function v = get.SomeProp(this)
            v = this.SomeProp;
        end
        
        function v = get.SomeDepProp(this)
            v = this.SomeProp * 5;
        end
    end
    
    events
        Test;
    end
    
end

