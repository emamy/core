classdef ICoreFun < dscomponents.IProjectable
    %ICOREFUN Basic interface for all dynamical system's core functions
    %
    % @Daniel Wirtz, 17.03.2010
    
    methods(Abstract)
        % Evaluates the core function
        y = evaluate(x,t,mu);
    end
    
end

