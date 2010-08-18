classdef classes_docGuide < handle
    %classes_docGuide Class summary - one line!
    %
    % This is the documentation guideline file for class documentation.
    %
    % < testing.MUnit & handle
    
    properties(SetAccess=private)
       
       % Some comment on the property.
       %
       % @type testing.MUnit
       %
       % See also: testing.MUnit
       SomeClass = testint.MUnit;
    end
    
    properties
       % Summary comment for SomeProp
       %
       % Detailed comment for SomeProp. Here you can write more detailed
       % text for the SomeProp property.
       %
       % @type int
       SomeProp = 0;
    end
    
    properties(Dependent)
       % equals someprop times five. 
       SomeDepProp;
    end
    

    methods
        function set.SomeProp(this, value)
            % Some function
            if value > 0
                this.SomeProp = value;
            else
                error('buh');
            end
        end
        
        
        function v = get.SomeProp(this)
            % Some redundant method's short description
            %
            % Some more comment!
            v = this.SomeProp;
        end
        
        function v = get.SomeDepProp(this)
            v = this.SomeProp * 5;
        end
    end
    
    methods(Sealed,Access=private)
        function noRealArguments(this)
            % This is the function brief description.
            %
            % Here are more details on the no real arguments function.
            % And even some more!
        end
    end
    
%     events
%         % This is the Test event's commentary.
%         Test;
%     end
    
end

