classdef ACoreFun < dscomponents.IProjectable
    % Basic interface for all dynamical system's core functions
    % Inherits the IProjectable interface.
    %
    % Subclassers have to implement the clone method according to the rules
    % explained in ICloneable, but only must re-implmement the project
    % method if a custom projection process takes place, i.e. coefficients
    % must be recomputed or such. Otherwise, the project implementation in
    % this method is sufficient to wrap 
    %
    % @author Daniel Wirtz @date 17.03.2010
    %
    % @change{0,3,sa,2011-04-15} Implemented Setters for the properties of this class 
    %
    % @change{0,3,dw,2011-04-13} the evaluate function supports multiargument evaluation now.
    %
    % @new{0,3,dw,2011-04-12} Added a new set-protected property
    % dscomponents.ACoreFun.MultiArgumentEvaluations. Allows for speedup in case the target function
    % can be called with a matrix of vectors instead of a single vector only.
    %
    % @todo Add test for non-custom projection (model with `W'f(x,t,\mu)V`)
    
    properties(Access=protected)
        % Set this property if the projection process is customized by
        % overriding the default project method. This affects the
        % evaluation method of the ACoreFun.
        CustomProjection = false;
    end
    
    properties(SetAccess=protected)
        % Flag that indicates whether this core function can be evaluated passing a matrix of
        % vectors instead of a single vector
        %
        % In case the custom implementation allows for such a call a significant speedup can be
        % achieved when using matlabs vector-based expressions
        %
        % @note Regarding parameters and their processing in multiargument-evaluations, the bsxfun
        % method from matlab may be of great help. See the models.pcd.CoreFun1D class for an
        % example.
        %
        % @default false
        %
        % See also: evaluate bsxfun models.pcd.CoreFun1D
        MultiArgumentEvaluations = false;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        V;
        W;
    end

    methods
        
        function target = project(this, V, W, target)
            % Sets the protected V,W matrices that can be utilized on core
            % function evaluations after projection.
            %
            % If the customized subclass sets the CustomProjection property
            % to true, this method must be called with a 3rd argument, the
            % target instance of a subclass upon which the projection
            % process is to be performed.
            %
            % Parameters:
            % V: Projection matrix one
            % W: Projection matrix two
            % target: The target instance which is a clone of this current
            % class. This parameter is required if custom projection is
            % used.
            %
            % See also: CustomProjection
            %
            % @change{0,3,dw,2011-04-11} Modified the project method to
            % also cater for the noncustomized projection case by cloning
            % the current instance (dscomponents.IProjectable now inherits
            % from ICloneable per default).
            if nargin < 4
                if this.CustomProjection
                    error('The target parameter must be given if custom projection is used.');
                end
                target = this.clone;
            end
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
            % mu: The parameter(s) to use. Set to [] if the function does not
            % support parameters.
            proj = ~(this.CustomProjection || isempty(this.V) || isempty(this.W));
            if proj
                x = this.V*x;
            end
            % check if fast evaluation is possible
            if size(x,2) == 1 || this.MultiArgumentEvaluations
                fx = this.evaluateCoreFun(x, t, mu);
            else
                % evaluate each point extra
                fx = zeros(size(x));
                for idx = 1:size(x,2)
                    fx(:,idx) = this.evaluateCoreFun(x(:,idx), t(idx), mu(:,idx));
                end
            end
            if proj
                fx = this.W'*fx;
            end
        end

        function set.CustomProjection(this, value)
            if ~islogical(value)
                error('Property must be logical/boolean. either true or false');
            end
            this.CustomProjection = value;
        end
        
        function set.MultiArgumentEvaluations(this, value)
            if ~islogical(value)
                error('Property must be logical/boolean. either true or false');
            end
            this.MultiArgumentEvaluations = value;
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

