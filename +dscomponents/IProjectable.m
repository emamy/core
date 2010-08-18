classdef IProjectable < handle
    % IPROJECTABLE Interface for all components that can be projected.
    %
    %  Projection in the KerMor context means restriction of the system's
    %  state space to a lower dimension. Therefore an orthogonal Matrix `V`
    %  is used and the system's components are multiplied from the left,
    %  right or even both sides by `V` or `V^t` respectively.
    %  However, for some components the projection process infers a change
    %  to the components' underlying data. Thus, when calling project a NEW
    %  instance of the implementing class has to be returned for which all
    %  functions behave according to the changes caused by the projection.
    %
    %  A special case is the approx.BaseApprox class which also implements
    %  this interface. Since the approx classes approximate the system's
    %  core function they will have approximation specific data that may
    %  change during the projection. Especially here a call to project must
    %  return a new instance with the projected approximation data (if the
    %  approximation method allows for it, of course).
    %
    % See also: approx BaseApprox
    %
    % @author Daniel Wirtz @date 17.03.2010
    % @docupdate
        
    methods(Abstract)
        % Returns a NEW INSTANCE of the projected object that does not rely
        % on data of the old one via references (everything must be copied
        % to ensure separability of reduced(=projected) versions and full
        % versions.
        % Parameters:
        % target: Optional (the target instance)
        %@docupdate
        projected = project(this, V, W, target);
    end
    
end