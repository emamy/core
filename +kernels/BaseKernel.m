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
        
        function fcn = getLipschitzFunction(this)
            % Method that allows error estimators to obtain a lipschitz
            % constant estimation function from this kernel.
            % 
            % The default is simply each kernel's global lipschitz constant
            % function. However, subclasses may override this method in
            % order to return a better (maybe local) lipschitz constant
            % estimation function. See the BellFunction implementation, for
            % example.
            %
            % See also: kernels.BellFunction error.BaseEstimator
            fcn = @this.getGlobalLipschitz;
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
        
        % Computes the partial derivatives with respect to each component of the first argument.
        %
        % Parameters:
        % x: The point where to evaluate the partial derivatives. Must be a single column `d\times 1` vector.
        % y: The corresponding center points at which the partial derivatives with respect to the
        % first argument are to be computed. Can be either a column vector `d\times 1` or a matrix `d\times n` containing
        % `n` multiple centers.
        %
        % Return values:
        % Nabla: A `d \times n` matrix of partial derivatives with respect to the first argument
        % evaluated using all second arguments.
        Nabla = getNabla(x,y)
        
        c = getGlobalLipschitz(this);
    end
    
end

