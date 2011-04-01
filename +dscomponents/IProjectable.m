classdef IProjectable < handle
    % IPROJECTABLE Interface for all components that can be projected.
    %
    %  Projection in the KerMor context means restriction of the system's
    %  state space to a lower dimension. Therefore an biorthogonal Matrix
    %  pair `V,W` is used and the system's components are multiplied from
    %  the left, right or even both sides by `V` or `W^t` respectively.
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
    %
    % @change{0,3,dw,2011-04-01} Updated documentation.
        
    methods(Abstract)
        % Returns a NEW INSTANCE of the projected object that does not rely
        % on data of the old one via references (everything must be copied
        % to ensure separability of reduced(=projected) versions and full
        % versions.
        %
        % Parameters:
        % V: The `V` matrix of the biorthogonal pair `V,W`
        % W: The `W` matrix of the biorthogonal pair `V,W`
        % target: Optional (the target instance)
        %
        % Return values:
        % projected: A new cloned instance with projection performed.
        %
        projected = project(this, V, W, target);
    end
    
end