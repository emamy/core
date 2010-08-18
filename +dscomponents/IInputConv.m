classdef IInputConv < dscomponents.IProjectable
    %BASEINPUTCONV Base class for input conversion "B".
    %   For simpler input conversions, it will be convenient to simply use
    %   the Pointer versions and pass the target function. For more complex
    %   input calculations which require local setup for example subclass
    %   this class and implement the evaluate method.
    %
    % See also: ACoreFun IOutputConv
    %
    % @author Daniel Wirtz @date 17.03.2010
    
    methods(Abstract)
        C = evaluate(t,mu);
    end
    
end

