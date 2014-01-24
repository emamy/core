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
        function this = PointerCoreFun(funPtr, xdim, multieval, timedep)
            % Creates a new core function `\vf` using the function pointer as inner function.
            %
            % Parameters:
            % funPtr: The function handle to use. Must implement the interface `(x,t,\mu)` to
            % be called with. @type function_handle
            % xdim: The input space dimension of the function `\vf`
            % multieval: If the function handle can take matrix valued `(x,t,mu)`, set this
            % flag to true to speed up computations. @type logical @default false
            % timedep: A flag that indicates if the passed function handle is (directly)
            % time-dependent or not. @type logical @default true
            %
            % @change{0,6,dw,2012-01-19} Added a new optional \c multieval parameter to
            % indicate that the function handle can take matrix valued arguments `x,t,\mu`.
            if nargin < 4
                % Assume worst case: set time dependency to true!
                timedep = true;
                if nargin < 3
                    multieval = false;
                end
            end
            % Creates a new wrapper for a core function handle.
            if ~isa(funPtr,'function_handle')
                error('Argument funPtr must be a function handle.');
            elseif nargin(funPtr) ~= 3
                error('funPtr nargin must equal three (= x,t,mu).');
            end
            this.CustomProjection = false;
            this.MultiArgumentEvaluations = multieval;
            this.TimeDependent = timedep;
            this.target = funPtr;
            this.xDim = xdim;
            this.fDim = size(funPtr(rand(xdim,1),0,rand(100,1)),1);
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
            proj = dscomponents.PointerCoreFun(newfun, size(V,2));
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
            m.System.x0 = dscomponents.ConstInitialValue(1);
            m.System.f = dscomponents.PointerCoreFun(fun,1);
            m.simulate();
            clear m;
        end
    end
    
end

