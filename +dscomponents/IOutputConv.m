classdef IOutputConv < dscomponents.IProjectable
    %BASEOUTPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
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

