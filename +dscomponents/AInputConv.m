classdef AInputConv < KerMorObject & general.AProjectable
    %AInputConv: Base class for input conversion "B".
    %   For simpler input conversions, it will be convenient to simply use
    %   the Pointer versions and pass the target function. For more complex
    %   input calculations which require local setup for example subclass
    %   this class and implement the evaluate method.
    %
    % See also: ACoreFun AOutputConv
    %
    % @author Daniel Wirtz @date 17.03.2010
    %
    % @change{0,5,dw,2011-07-07} Fixed output name from `C` to `B`
    
    properties(SetObservable)
        % Flag that indicates if the AInputConv is (truly) time-dependent.
        %
        % Set in subclasses to the correct value for your implementation.
        % Most built-in KerMor classes set their values correctly or guess them.
        %
        % This will cause different behaviour for e.g. ODE solvers.
        %
        % @propclass{critical} Not setting this value in implementing subclasses causes
        % KerMor's ODE solvers to (possibly) produce wrong results due to wrong assumptions on
        % the time dependence of the core function.
        %
        % @type logical @default true
        TimeDependent = true;
    end
    
    methods
        function prepareSimulation(this, mu)
            % do nothing by default
        end
    end
    
    methods(Abstract)
        % Template method that evaluates the input conversion matrix `B` at the current time `t`
        % and [optional] parameter `\mu`.
        B = evaluate(t, mu);
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, varargin)
            obj = loadobj@KerMorObject(obj, varargin{:});
            obj = loadobj@general.AProjectable(obj, varargin{:});
        end
    end
    
end

