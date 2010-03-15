classdef ICoreFun < dscomponents.IProjectable
    %BASECOREFUN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Abstract)
        y = evaluate(x,t,mu);
    end
    
end

