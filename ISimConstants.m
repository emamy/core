classdef ISimConstants < handle
    %ISIMCONSTANTS Interface for systems, core functions or other components that have
    % properties that are constant during a simulation.
    %
    % See also: BaseDynSystem ACoreFun
    %
    % @author Daniel Wirtz @date 2010-03-17
    %
    % @change{0,4,dw,2011-05-31} Generalized this interface description so it could be used for any
    % component. (New: error.BaseLocalLipschitzFunction)
    
    methods(Abstract)
        % Initializes inner properties that stay constant for the
        % simulation of a trajectory. 
        %
        % This is a chance to gain efficiency for some cases, for example
        % space-discretization variables that are to be computed for a
        % simulation but not every timestep.
        %
        % See also: BaseDynSystem ACoreFun BaseLocalLipschitzFunction
        prepareConstants;
    end
    
end

