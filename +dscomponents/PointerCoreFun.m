classdef PointerCoreFun < dscomponents.ICoreFun
    %POINTERCOREFUN Allows for core functions provided by function handles.
    %   
    % In many contexts the creation of a specific class implementing
    % ICoreFun is not necessary due to the simplicity of the core function
    % or because no function-specific properties have to be changed or
    % modeled. In this case use this class and pass a function handle to
    % the constructor which will be used as the actual core function.
    %
    % See also: PointerInputConv PointerOutputConv
    %
    % @Daniel Wirtz, 16.03.2010
    
    properties(Access=private)
        Target;
    end
    
    methods
        function this = PointerCoreFun(funPtr)
            % Creates a new wrapper for a core function handle.
            if ~isa(funPtr,'function_handle')
                error('Argument funPtr must be a function handle.');
            elseif nargin(funPtr) ~= 3
                error('funPtr nargin must equal three (= x,t,mu).');
            end
            this.Target = funPtr;
        end
        
        function updateSimConstants(this)%#ok
            % Nothing to do here.
        end
        
        function fx = evaluate(this, x, t, mu)
            % Evaluates the core function at x,t,mu
            fx = this.Target(x,t,mu);
        end
        
        function proj = project(this, V)
            % Projects the core function into the reduced space.
            % Creates a new PointerCoreFun and computes `\hat{f}(z) =
            % V^tf(Vz)`
            newfun = @(z,t,mu)V' * this.Target(V*z,t,mu);
            proj = dscomponents.PointerCoreFun(newfun);
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

