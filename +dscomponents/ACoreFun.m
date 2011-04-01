classdef ACoreFun < dscomponents.IProjectable & ICloneable
    % Basic interface for all dynamical system's core functions
    % Inherits the IProjectable interface.
    %
    % @author Daniel Wirtz @date 17.03.2010
    %
    % @todo Add test for non-custom projection (model with `W'f(x,t,\mu)V`)
    
    properties(Access=protected)
        % Set this property if the projection process is customized by
        % overriding the default project method. This affects the
        % evaluation method of the ACoreFun.
        CustomProjection = false;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        V;
        W;
    end

    methods
        
        function target = project(this, V, W, target)%#ok
            % Sets the protected V,W matrices that can be utilized on core
            % function evaluations after projection.
%             if nargin == 4
%                 projected = target;
%             else
%                 projected = this.clone;
%             end
%             projected.V = V;
%             projected.W = W;
            target.V = V;
            target.W = W;
        end
        
        function fx = evaluate(this, x, t, mu)
            % Evaluates the f-approximation. Depending on a possible
            % projection and the CustomProjection-property the function
            % either calls the inner evaluation directly which assumes
            % `f = f^r(z)` or projects the reduced state variable z into
            % the original space and evaluates the function there, so via
            % `f = V'f(Vz)`
            %
            % Parameters:
            % x: The state variable vector/matrix (with colum state
            % vectors)
            % t: The corresponding times for each state vector. Set to []
            % if no time is used.
            % mu: The parameter to use. Set to [] if the function does not
            % support parameters.
            if this.CustomProjection || isempty(this.V) || isempty(this.W)
                fx = this.evaluateCoreFun(x, t, mu);
            else
                fx = this.W'*this.evaluateCoreFun(this.V*x, t, mu);
            end
        end

        function set.CustomProjection(this, value)
            if ~islogical(value)
                error('Property "CustomProjection" must be logical/boolean.');
            end
            this.CustomProjection = value;
        end
        
        function copy = clone(this, copy)
            if nargin == 1 || ~isa(copy,'dscomponents.ACoreFun')
                error('Incorrect call to clone. As this class is abstract, a subclass of ACoreFun has to be passed as second argument.');
            end
            % Copy local properties
            copy.CustomProjection = this.CustomProjection;
            copy.V = this.V;
            copy.W = this.W;
        end
    end
        
    methods(Abstract)
        % Evaluates the core function
        y = evaluateCoreFun(this, x, t, mu);
    end
    
end

