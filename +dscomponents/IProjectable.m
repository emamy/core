classdef IProjectable < handle
    
    methods(Abstract)
        projected = project(V);
    end
    
end