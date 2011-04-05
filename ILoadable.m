classdef ILoadable < handle
    % Interface for any objects that implement customized loading
    % behaviour via MatLab's built-in loadobj-method.
    %
    % Type @code >> docsearch save and load process @endcode at the matlab
    % console for more information or visit
    % http://www.mathworks.com/help/techdoc/matlab_oop/bres1y6.html
    % to get information about how matlab interacts with those methods.
    %
    % Loading guidelines:
    % -# Function signature
    %  -# The loadable class may be subclassed later: 
    %    The method must take an extra argument "obj" which contains a
    %    subclass instance. Then all local attributes (public
    %    and private) must be set at that passed instance. This
    %    ensures that all properties along the inheritance tree are
    %    taken care of.
    %  -# The class is a final class / last in the inheritance chain: 
    %    No extra argument is necessary, just create a new blank
    %    instance 'obj' using the classes constructor (including any
    %    necessary constructor arguments; these will be derivable from the
    %    classes context)
    % -# Call to superclass clone method
    %   For subclasses you have to call the loadobj method of the next
    %   superior class(es) via @code obj = loadobj@superclassname(s, obj)
    %   @endcode
    % -# Setting properties
    %   Only assign properties that are declared within that very same
    %   class. Take care that you assign those properties in an order that
    %   avoid some setting methods to throw errors as necessary properties
    %   have not been set yet.
    %
    % See also: load save ICloneable
    %
    % @author Daniel Wirtz @date 2011-04-05
    %
    % @new{0,3,dw,2011-04-05} Added this class to support KerMor-wide
    % correct loading of objects and models.
    
    methods(Static, Access=protected)
        % Performs custom object loading with user-defined order of
        % property assignment.
        %
        % See the classes documentation for more information about how to
        % work with this loadobj method.
        %  
        % This method is called my MatLab automatically upon loading an
        % instance of this class or subclasses, unless it is overridden
        % in any subclass. The last is necessary for any subclasses to
        % implement in order to successfully load saved instances,
        % since at this point the actual subclass and its properties
        % cannot be created properly (this is to be done inside a
        % loadobj-method override inside the subclass)
        %
        % Any classes in KerMor that are abstract and implement a
        % loadobj-method check for the second argument to be passed. If
        % one subclasses one of those classes and does not implement a
        % loadobj which calls (all) the superclasses' loadobj methods
        % passing the concrete instance, an error will be issued.
        %
        % The reason why every abstract class checks for that
        % additional parameter is simply that those classes can be
        % subclassed at any point in the hierarchy, depending on the
        % customers needs.
        obj = loadobj(s, obj);
    end
    
end

