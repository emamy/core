classdef PointerCoreFun < dscomponents.ACoreFun
    % Allows for core functions provided by function handles.
    %   
    % In many contexts the creation of a specific class implementing
    % ACoreFun is not necessary due to the simplicity of the core function
    % or because no function-specific properties have to be changed or
    % modeled. In this case use this class and pass a function handle to
    % the constructor which will be used as the actual core function.
    %
    % See also: PointerInputConv PointerOutputConv
    %
    % @author Daniel Wirtz @date 16.03.2010
    
    properties(Access=private)
        target;
    end
    
    methods
        function this = PointerCoreFun(funPtr, multieval)
            % Creates a new core function using the function pointer as inner function.
            %
            % Parameters:
            % funPtr: The function handle to use. Must implement the interface `(x,t,\mu)` to
            % be called with. @type function_handle
            % multieval: If the function handle can take matrix valued `(x,t,mu)`, set this
            % flag to true to speed up computations. @type logical @default false
            %
            % @change{0,6,dw,2012-01-19} Added a new optional \c multieval parameter to
            % indicate that the function handle can take matrix valued arguments `x,t,\mu`.
            if nargin == 1
                multieval = false;
            end
            % Creates a new wrapper for a core function handle.
            if ~isa(funPtr,'function_handle')
                error('Argument funPtr must be a function handle.');
            elseif nargin(funPtr) ~= 3
                error('funPtr nargin must equal three (= x,t,mu).');
            end
            this.CustomProjection = false;
            this.MultiArgumentEvaluations = multieval;
            this.target = funPtr;
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            % Evaluates the core function at x,t,mu
            fx = this.target(x,t,mu);
        end
                
        function proj = project(this, V, W)
            % Projects the core function into the reduced space.
            % Creates a new PointerCoreFun and computes `\hat{f}(z) =
            % V^tf(Vz)`
            
            % no need to call clone here as the only property gets changed
            % anyways.
            newfun = @(z,t,mu)W' * this.target(V*z,t,mu);
            proj = dscomponents.PointerCoreFun(newfun);
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = dscomponents.PointerCoreFun(this.target);
        end
    end
    
    methods(Static)
        function test_PointerCoreFun
            m = models.BaseFullModel;
            fun = @(x,t,mu)x*3.*t;
            m.System.f = dscomponents.PointerCoreFun(fun);
            m.simulate();
            clear m;
        end
    end
    
end

