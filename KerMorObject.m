classdef KerMorObject < handle
    % Base class for any KerMor class
    %
    % Any KerMor object implements custom loading behaviour.
    %
    % Moreover, a messaging system is about to be implemented in KerMor
    % to obtain different modes of output and logging. Every class should
    % have access to that system, therefore a protected method will be
    % available within any subclass. (Current skeleton: logMessage)
    % 
    % @author Daniel Wirtz @date 2011-04-05
    %
    % @change{0,3,dw,2011-04-07} 
    % - Moved the checkType method from models.BaseModel to this class.
    % - Removed the ALoadable class and moved the ID property here.
    %
    % @new{0,3,dw,2011-04-05} Added this class to give any object a common
    % base class and some shared functionality. For this version it will be
    % loading and message logging, possibly cloning as well.
    
    properties(SetAccess=private)
        % An ID that allows to uniquely identify this KerMorObject
        ID;
    end
    
    methods
        function this = KerMorObject
            this.ID = general.IDGenerator.generateID;
        end
    end
    
    methods(Access=protected)
        function logMessage(this, msg)
            % Puts a message to the KerMor messaging & logging system
            
%             % maybe: notify any event listeners about the logged message?
%             data = event.EventData;
%             notify(this,'Log',data);
        end
    end
    
    methods(Access=protected,Sealed)
        function checkType(this, obj, type)%#ok
            % Object typechecker.
            %
            % Checks if a given object is of the specified type and throws
            % an error if not. Convenience method.
            if ~isempty(obj) && ~isa(obj, type)
                error(['Wrong type ''' class(obj) ''' for this property. Has to be a ' type]);
            end
        end
    end

%     events
%         Log;
%     end
    
end

