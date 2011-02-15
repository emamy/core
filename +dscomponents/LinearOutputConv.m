classdef LinearOutputConv < dscomponents.AOutputConv
    %STDOUTPUTCONV Standard linear output converter.
    %
    %   Simply forwards the system's states to the output.
    %   If projected, the states are projected back into the original sized
    %   state space.
    %
    % @author Daniel Wirtz @date 15.03.2010
    
    properties(Access=private)
        C;
    end
    
    methods
        function this = LinearOutputConv(C)
            % Standard constructor.
            this.C = C;
            
            % Ensure that the flag is set correctly (independent from the
            % default setting which may change)
            this.TimeDependent = false;
        end
        
        function copy = project(this, V, W)%#ok
            % Performs projection for the standard output conversion.
            %
            % Simply creates a new instance (to keep the full model
            % running) and sets the projection matrix.
            copy = dscomponents.LinearOutputConv(this.C*V);
        end
        
        function C = evaluate(this, t, mu)%#ok
            % Evaluates the output conversion matrix.
            % In this simple case this is just the projection matrix, if
            % set, otherwise 1.
            C = this.C;
        end
    end
    
end

