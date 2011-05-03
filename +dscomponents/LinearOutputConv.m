classdef LinearOutputConv < dscomponents.AOutputConv
    % Standard linear output converter.
    %
    % Simply forwards the system's states to the output.
    % If projected, the states are projected back into the original sized
    % state space.
    %
    % @author Daniel Wirtz @date 15.03.2010
    %
    % @change{0,3,dw,2011-04-11} Due to the new inheritance IProjectable <
    % ICloneable the clone method is implemented here now, too.
    
    properties(Access=private)
        C;
    end
    
    methods
        function this = LinearOutputConv(C)
            this = this@dscomponents.AOutputConv;
            
            % Standard constructor.
            this.C = C;
            
            % Ensure that the flag is set correctly (independent from the
            % default setting which may change)
            this.TimeDependent = false;
        end
        
        function copy = project(this, V, W)%#ok
            % Performs projection for the standard output conversion.
            copy = this.clone;
            copy.C = this.C*V;
        end
        
        function copy = clone(this)
            copy = dscomponents.LinearOutputConv(this.C);
            copy = clone@dscomponents.AOutputConv(this, copy);
        end
        
        function C = evaluate(this, t, mu)%#ok
            % Evaluates the output conversion matrix.
            % In this simple case this is just the projection matrix, if
            % set, otherwise 1.
            C = this.C;
        end
    end
    
end

