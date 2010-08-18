classdef BaseKernel < handle
    %BASEKERNEL Basic KerMor Kernel
    %
    % All Kernels have to inherit from this class.
    % @todo remove RotInv property and replace for check of interface
    % IRotInv implementation instead
    %
    % @author Daniel Wirtz @date 12.03.2010
    
    properties(SetAccess = protected)
        % Flag whether this kernel is rotation-invariant.
        RotationInvariant;
    end
    
    methods
        function this = BaseKernel
            this.RotationInvariant = isa(this, 'kernels.IRotationInvariant');
        end
    end
        
    methods(Abstract)
        % Evaluation method for the current kernel.
        %
        % Convention:
        % The second parameter is OPTIONAL, if not passed x=y is assumed.
        % For reasons of possible efficiency increase we let the
        % implementer handle the check.
        K = evaluate(x,y);
        
        c = getGlobalLipschitz(this);
    end
    
end

