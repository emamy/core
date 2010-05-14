classdef IOutputConv < dscomponents.IProjectable
    %BASEOUTPUTCONV Base class for output conversion "C".
    %   For simpler output conversions, it will be convenient to simply use
    %   the Pointer versions and pass the target function. For more complex
    %   output calculations which require local setup for example subclass
    %   this class and implement the evaluate method.
    %
    % See also: ICoreFun IInputConv
    %
    % @DanielWirtz, 17.03.2010
    
    properties(SetAccess=protected)
        % Flag whether the output converter actually depends on a time
        % variable. 
        % Implemented for speed reasons when computing the output. This
        % flag can (and should!) be set in classes that implement this
        % interface.
        %
        % Defaults to false.
        TimeDependent = false;
    end
    
    methods(Abstract)
        y = evaluate(t,mu);
    end
    
end

