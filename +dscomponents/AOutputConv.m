classdef AOutputConv < dscomponents.IProjectable
    %BASEOUTPUTCONV Base class for output conversion "C".
    %   For simpler output conversions, it will be convenient to simply use
    %   the Pointer versions and pass the target function. For more complex
    %   output calculations which require local setup for example subclass
    %   this class and implement the evaluate method.
    %
    % See also: ACoreFun IInputConv
    %
    % @author Daniel Wirtz @date 17.03.2010
    
    properties(SetAccess=protected)
        % Flag whether the output converter actually depends on a time
        % variable. 
        % Implemented for speed reasons when computing the output. This
        % flag can (and should!) be set in classes that implement this
        % interface.
        %
        % Defaults to false.
        TimeDependent = false;
    end
    
    methods
        function y = computeOutput(this, t, x, mu)
            % Computes the output `y(t) = C(t,\mu)x(t)` from a given state
            % result vector `x(t)`.
            %
            % Convenience method, as it is currently used in
            % models.BaseModel/simulate and error.DefaultEstimator/process
            %
            % Parameters:
            % t: The time-step vector
            % x: The state variable vector at each time step per column
            % mu: The current parameter mu. Optional.
            %
            % Return values:
            % y: The output according to `y(t) = C(t,\mu)x(t)`
            if nargin == 2
                mu = [];
            end
            if this.TimeDependent
                % Evaluate the output conversion at each time t
                % Figure out resulting size of C*x evaluation
                hlp = this.evaluate(t(1),mu)*x(:,1);
                y = zeros(size(hlp,1),length(t));
                y(:,1) = hlp;
                for idx=2:length(t)
                    y(:,idx) = this.evaluate(t(idx),mu)*x(:,idx);
                end
            else
                % otherwise it's a constant matrix so multiplication
                % can be preformed much faster.
                y = this.evaluate([],mu)*x;
            end
        end
        
        function copy = clone(this, copy)
            copy.TimeDependent = this.TimeDependent;
        end
    end
    
    methods(Abstract)
        y = evaluate(t,mu);
    end
    
end