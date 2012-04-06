classdef BaseKernel < KerMorObject
    %BASEKERNEL Basic KerMor Kernel
    %
    % All Kernels have to inherit from this class.
    %
    % @author Daniel Wirtz @date 12.03.2010
    %
    % @change{0,3,dw,2011-04-21} Removed the RotationInvariant property as it is now replaced by the
    % IRotationInvariant interface.
    %
    % @todo 
    % - Implement cloning for kernels
    % - write universal tests for kernels that check the interface functions (IJacobian etc)
    % - Implement eq-method for all kernels
    
    properties(SetObservable)
        % Projection/selection matrix for argument components
        %
        % Set this value to the indices of the components of any argument passed to the kernel
        % that should be effectively used. This property is mainly used with parameter kernels
        % in order to extract relevant entries. Leave to [] if all values should be used.
        %
        % Subclasses must take care to use this property if set.
        %
        % @propclass{data} Depends on the kernel setting and problem setup.
        %
        % @default [] @type matrix
        P = [];
    end
    
    methods
        function this = BaseKernel
            this = this@KerMorObject;
            this.registerProps('P');
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
        
        function bool = eq(A, B)
            % Checks if a kernel equals another kernel
            bool = eq@KerMorObject(A,B) && isequal(A.P, B.P);
        end
    end
        
    methods(Abstract)
        % Evaluation method for the current kernel.
        %
        % Parameters:
        % x: First set `x_i \in \R^d` of `n` vectors @type matrix<double>
        % y: Second set `y_j \in \R^d` of `m` vectors. If y is empty `y_i = x_i` and `n=m`
        % is to be assumed. @type matrix<double>
        %
        % Return values:
        % Phi: The evaluation matrix `\Phi(x,y) \in \R^{n\times m}` of the kernel `\Phi`, with
        % entries `\Phi(x_i,y_j)` at `i,j`.
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
        Nabla = getNabla(this, x, y)
        
        % Returns the global lipschitz constant of this kernel.
        %
        % Exprimental state as not implemented & checked for all kernels.
        c = getGlobalLipschitz(this);
    end
    
end

