classdef AInputConv < KerMorObject & dscomponents.IProjectable
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
    
    methods(Abstract)
        % Template method that evaluates the input conversion matrix `B` at the current time `t`
        % and [optional] parameter `\mu`.
        B = evaluate(t, mu);
    end
    
end

