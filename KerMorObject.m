classdef KerMorObject < ILoadable
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
    % @new{0,3,dw,2011-04-05} Added this class to give any object a common
    % base class and some shared functionality. For this version it will be
    % loading and message logging, possibly cloning as well.
    
    methods(Access=protected)
        function logMessage(this, msg)
            % Puts a message to the KerMor messaging & logging system
            
%             % maybe: notify any event listeners about the logged message?
%             data = event.EventData;
%             notify(this,'Log',data);
        end
    end
    
%     events
%         Log;
%     end
    
end

