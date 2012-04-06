classdef KerMorObject < DPCMObject
    % Base class for any KerMor class
    %
    % Every KerMor object is identifiable by a unique ID.
    %
    % Moreover, a messaging system is about to be implemented in KerMor
    % to obtain different modes of output and logging. Every class should
    % have access to that system, therefore a protected method will be
    % available within any subclass. (Current skeleton: logMessage)
    % 
    % @author Daniel Wirtz @date 2011-04-05
    %
    % @new{0,3,dw,2011-04-21} Implemented the property default changed supervision system as
    % described in @ref propclasses. Those changes affect most KerMor classes.
    %
    % @change{0,3,dw,2011-04-07} 
    % - Moved the checkType method from models.BaseModel to this class.
    % - Removed the ALoadable class and moved the ID property here.
    %
    % @new{0,3,dw,2011-04-05} Added this class to give any object a common
    % base class and some shared functionality. For this version it will be
    % loading and message logging, possibly cloning as well.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    %
    % @todo Extract DPCS into single class and put onto matlab file exchange!
    
    methods
        function this = KerMorObject
            % Constructs a new KerMor object.
            %
            % Important notice at this stage: Due to possible multiple inheritance any object's
            % constructor should check if any custom properties that are assigned during
            % construction are already present / different from the default value given at property
            % declaration. If so, chances are that the constructor is called again for the same
            % object; in this case, assigning a new ID and property changed dictionary caused any
            % old registered properties to be overwritten!
            
            this = this@DPCMObject;
        end
        
        function display(this)
            %disp(object2str(this));
            disp(this);
        end
        
        function bool = eq(A, B)
            % Checks equality of two KerMor objects.
            %
            % Override (with superclass method call) in subclasses for
            % specific comparison.
            bool = strcmp(class(A), class(B));
        end
        
        function bool = ne(A, B)
            % Checks if two KerMorObjects are different
            %
            % Negates the result of eq(A,B) for KerMorObject.
            bool = ~eq(A,B);
        end
    end
        
    methods(Sealed, Access=protected)
        function checkType(this, obj, type)%#ok
            % Object typechecker.
            %
            % Checks if a given object is of the specified type and throws
            % an error if not. Convenience method that can bu used in subclasses at any setter
            % method that ensures the correct value type.
            if ~isempty(obj) && ~isa(obj, type)
                error(['Wrong type ''' class(obj) ''' for this property. Has to be a ' type]);
            end
        end
    end
end

