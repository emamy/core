classdef BaseSpaceReducer < handle
    % Base class for all space reduction algorithms.
    %
    % @author Daniel Wirtz
    % @date 11.03.2010
    
    properties
        % Maybe some prop like ndims here?
    end
    
    methods(Abstract)
        [V,W] = generateReducedSpace(model);
    end
    
end

