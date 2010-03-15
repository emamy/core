classdef ICloneable < handle
    %ICLONEABLE Interface for cloneable handle classes
    %
    % Any implementing class should optionally take the "target" argument
    % as an subclass instance of the current class. This is to enable
    % superclasses to fill in private variables needed for correct
    % functionality of the cloned object.
    % Example: Superclass "A" has a private property "d", and subclass "B"
    % is calling A.clone(<new B instance> new) so that class "A" can set
    % new.d = this.d.
    %
    % Cloning guidelines:
    % -# Function signature
    %  - The cloneable class may be subclassed:
    %    The method must take an extra argument "target" which may
    %    contain a subclass instance. Then all local attributes (public
    %    and private) must be copied into that passed instance. This
    %    ensures that all properties along the inheritance tree are
    %    taken care of.
    %  - The class is a final class:
    %    No extra argument is necessary, just create a new blank
    %    instance using the class' constructor (including any necessary
    %    constructor arguments; these will be derivable from the
    %    classes' context)
    % -# Call to superclass clone method 
    %   For subclasses you have to call the clone method of the next
    %   superior class via @code clone@superclassname(this, target)
    %   @endcode
    % -# Properties
    %   Only copy properties that are declared within that very same class!
    %
    % @Daniel Wirtz, 15.03.2010
    
    methods(Access=protected, Abstract)
        copy = clone(target);
    end
    
end

