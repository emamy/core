classdef PointerInputConv < dscomponents.AInputConv
    %POINTERINPUTCONV Allows for input converters provided by function handles.
    %   
    % In many contexts the creation of a specific class implementing
    % AInputConv is not necessary due to the simplicity of the conversion
    % or because no function-specific properties have to be modeled. In
    % this case use this class and pass a function handle to the
    % constructor which will be used as the actual function.
    %
    % See also: PointerCoreFun PointerOutputConv
    %
    % @author Daniel Wirtz @date 16.03.2010
    %
    % @change{0,3,dw,2011-04-11} Implemented the clone method from
    % ICloneable
    
    properties(Access=private)
        Target;
    end
    
    methods
        function this = PointerInputConv(funPtr)
            % Creates a new wrapper for a core function handle.
            if ~isa(funPtr,'function_handle')
                error('Argument funPtr must be a function handle.');
            elseif nargin(funPtr) ~= 2
                error('funPtr nargin must equal two (= t,mu).');
            end
            this.Target = funPtr;
        end
        
        function B = evaluate(this, t, mu)
            % Evaluates the core function at t,mu
            B = this.Target(t,mu);
        end
      
        function proj = project(this, V, W)
            % Projects the core function into the reduced space.
            % Creates a new PointerInputConv and computes `B^r(t,\mu) =
            % V^tB(t,\mu)`
            newfun = @(t,mu)W' * this.Target(t,mu);
            proj = dscomponents.PointerInputConv(newfun);
            proj = project@general.AProjectable(this, V, W, proj);
        end
        
        function copy = clone(this)
            copy = dscomponents.PointerInputConv(this.Target);
        end
    end    
end

