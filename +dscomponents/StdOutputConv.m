classdef StdOutputConv < dscomponents.IOutputConv
    %STDOUTPUTCONV Standard output converter.
    %
    %   Simply forwards the system's states to the output.
    %   If projected, the states are projected back into the original sized
    %   state space.
    %
    % @author Daniel Wirtz @date 15.03.2010
    
    properties(Access=private)
        V;
    end
    
    methods
        
        function this = StdOutputConv
            % Standard constructor.
            
            % Ensure that the flag is set correctly (independent from the
            % default setting which may change)
            this.TimeDependent = false;
        end
        
        function copy = project(this, V, W)%#ok
            % Performs projection for the standard output conversion.
            %
            % Simply creates a new instance (to keep the full model
            % running) and sets the projection matrix.
            copy = dscomponents.StdOutputConv;
            copy.V = V;
        end
        
        function C = evaluate(this, t, mu)%#ok
            % Evaluates the output conversion matrix.
            % In this simple case this is just the projection matrix, if
            % set, otherwise 1.
            if ~isempty(this.V)
                C = this.V;
            else
                C = 1;
            end
        end
    end
    
end

