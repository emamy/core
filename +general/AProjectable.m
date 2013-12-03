classdef AProjectable < ICloneable
    %  Interface for all components that can be projected.
    %
    %  Projection in the KerMor context means restriction of the system's
    %  state space to a lower dimension. Therefore an biorthogonal Matrix
    %  pair `V,W` is used and the system's components are multiplied from
    %  the left, right or even both sides by `V` or `W^t` respectively.
    %  However, for some components the projection process infers a change
    %  to the components' underlying data. Thus, when calling project a NEW
    %  instance of the implementing class has to be returned for which all
    %  functions behave according to the changes caused by the projection.
    %  Due to this every projectable class also implements ICloneable in
    %  order to allow to create instance copies upon which the projection
    %  process can be performed.
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
    % @author Daniel Wirtz @date 2010-03-17
    %
    % @change{0,6,dw,2012-06-06} Moved from dscomponents to general package
    %
    % @change{0,3,dw,2011-04-11} Inheriting from ICloneable now in order to
    % ensure cloning capabilities for any projectable class.
    %
    % @change{0,3,dw,2011-04-01} Updated documentation.
    %
    % @todo change AProjectable to IProjectable and remove V,W properties (not needed in all
    % components that are projectable)
    
    properties(SetAccess=protected)
        % The `V` matrix of the biorthogonal pair `V,W`
        V;
        
        % The `W` matrix of the biorthogonal pair `V,W`
        W;
    end
        
    methods
        function target = project(this, V, W, target)
            % Returns a NEW INSTANCE of the projected object that does not rely
            % on data of the old one via references (everything must be copied
            % to ensure separability of reduced(=projected) versions and full
            % versions, unless.
            %
            % Override in inheriting classes and pass the subclasses cloned
            % instance as fourth target parameter.
            %
            % Parameters:
            % V: The `V` matrix of the biorthogonal pair `V,W` @type
            % matrix<double>
            % W: The `W` matrix of the biorthogonal pair `V,W` @type
            % matrix<double>
            % target: Specify the subclasses projection target instance if
            % project is overridden in a subclass and the subclass has been
            % cloned. @type handle @default []
            %
            % Return values:
            % projected: A new cloned instance with projection performed.
            if nargin < 4
                target = this.clone;
            end
            target.V = V;
            target.W = W;
        end
        
        function copy = clone(this, copy)
            if nargin == 1
                copy = general.AProjectable;
            end
            copy.V = this.V;
            copy.W = this.W;
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if nargin == 2
                obj.V = from.V;
                obj.W = from.W;
            end
        end
    end
    
end