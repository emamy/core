classdef PointerOutputConv < dscomponents.IOutputConv
    %POINTEROUTPUTCONV Allows for input converters provided by function handles.
    %   
    % In many contexts the creation of a specific class implementing
    % IOutputConv is not necessary due to the simplicity of the conversion
    % or because no function-specific properties have to be modeled. In
    % this case use this class and pass a function handle to the
    % constructor which will be used as the actual function.
    %
    % See also: PointerCoreFun PointerInputConv
    %
    % @Daniel Wirtz, 16.03.2010
    
    properties(Access=private)
        Target;
    end
    
    methods
        function this = PointerOutputConv(funPtr, time_dependent)
            % Creates a new wrapper for a core function handle.
            %
            % Parameters:
            % funPtr: A pointer to the target function.
            % time_dependent: Optional parameter that determines if the
            % output conversion function is time-dependent.
            if ~isa(funPtr,'function_handle')
                error('Argument funPtr must be a function handle.');
            elseif nargin(funPtr) ~= 2
                error('funPtr nargin must equal two (= t,mu).');
            end
            this.Target = funPtr;
            % only set if given. Default value see IOutputConv
            if nargin == 2
                this.TimeDependent = time_dependent;
            end
        end
        
        function B = evaluate(this, t, mu)
            % Evaluates the core function at t,mu
            B = this.Target(t,mu);
        end
        
        function proj = project(this, V)
            % Projects the core function into the reduced space.
            % Creates a new PointerOutputConv and computes `C^r(t,\mu) =
            % C(t,\mu)V`
            newfun = @(t,mu)this.Target(t,mu) * V;
            proj = dscomponents.PointerOutputConv(newfun, this.TimeDependent);
        end
    end
    
end



