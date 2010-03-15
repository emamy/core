classdef IInputConv < dscomponents.IProjectable
    %BASEINPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Abstract)
        C = evaluate(t,mu);
    end
    
end

