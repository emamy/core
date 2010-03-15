classdef BaseKernel < handle
    %BASEKERNEL Basic KerMor Kernel
    %
    % All Kernels have to inherit from this class.
    %
    % @Daniel Wirtz, 12.03.2010
    
    properties(SetAccess=protected)
        % Flag whether this kernel is rotation-invariant.
        RotationInvariant = false;
    end
    
%     methods
%         function K = evaluateIdent(this, x)
%             K = this.evaluate(x,x);
%         end
%     end
    
    methods(Abstract)
        % Evaluation method for the current kernel.
        %
        % Convention:
        % The second parameter is OPTIONAL, if not passed x=y is assumed.
        % For reasons of possible efficiency increase we let the
        % implementer handle the check.
        K = evaluate(x,y);
    end
    
end

