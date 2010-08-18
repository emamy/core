classdef ISimConstants < handle
    %ISIMCONSTANTS Interface for systems or core functions that have
    %properties that are constant during simulation.
    %
    % See also: BaseDynSystem ACoreFun
    %
    % @author Daniel Wirtz @date 17.03.2010
    
    methods(Abstract)
        % Initializes inner properties that stay constant for the
        % simulation of a trajectory. 
        %
        % This is a chance to gain efficiency for some cases, for example
        % space-discretization variables that are to be computed for a
        % simulation but not every timestep.
        %
        % See also: BaseDynSystem ACoreFun
        updateSimConstants;
    end
    
end

