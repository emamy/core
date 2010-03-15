classdef BaseSpaceReducer < handle
    %BASERSPACE Base class for all space reduction algorithms.
    %
    % @DanielWirtz, 11.03.2010
    
    properties
        % Maybe some prop like ndims here?
    end
    
    methods(Abstract)
        V = generateReducedSpace(model);
    end
    
end

