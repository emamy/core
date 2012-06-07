classdef AOutputConv < KerMorObject & general.AProjectable
    %BASEOUTPUTCONV Base class for output conversion "C".
    %   For simpler output conversions, it will be convenient to simply use
    %   the Pointer versions and pass the target function. For more complex
    %   output calculations which require local setup for example subclass
    %   this class and implement the evaluate method.
    %
    % See also: ACoreFun AInputConv
    %
    % @author Daniel Wirtz @date 17.03.2010
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @change{0,3,sa,2011-04-15} Implemented Setter for the class property 
    
    properties(SetObservable,SetAccess=protected)
        % Flag whether the output converter actually depends on a time
        % variable. 
        % Implemented for speed reasons when computing the output. This
        % flag can (and should!) be set in classes that implement this
        % interface.
        %
        % @propclass{critical} Some output conversion matrices are time dependent. This property
        % must be set to the correct value in order for the output conversion to work correctly.
        %
        % @type logical
        %
        % Defaults to false.
        TimeDependent = false;
    end
    
    methods
        function this = AOutputConv
            this = this@KerMorObject;
            this.registerProps('TimeDependent');
        end
        
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
        
        function set.TimeDependent(this, value)
            if ~islogical(value)
                error('Property is logical. Must be set either true or false');
            end
            this.TimeDependent = value;
        end
        
        function copy = clone(this, copy)
            copy.TimeDependent = this.TimeDependent;
        end
    end
    
    methods(Abstract)
        % Template method that evaluates the output conversion matrix `C` at the current time `t`
        % and [optional] parameter `\mu`.
        y = evaluate(t,mu);
    end
    
end