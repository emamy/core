classdef BaseKernel < KerMorObject & ICloneable
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
        % The matrix `\vG` that induces the state space scalar product
        % `\spG{x}{y}` and norm `\noG{x-y}` to use.
        %
        % Must be a positive definite, symmetric matrix `\vG`.
        %
        % @propclass{critical} If a custom norm is used (i.e. after
        % subspace projection) tihs must be set in order to obtain correct
        % evaluations.
        %
        % @type matrix<double> @default 1
        G = 1;
        
        % Projection/selection matrix `\vP` for argument components
        %
        % Set this value to the indices of the components of any argument passed to the kernel
        % that should be effectively used. This property is mainly used with parameter kernels
        % in order to extract relevant entries. Leave to [] if all values should be used.
        %
        % Subclasses must take care to use this property if set.
        %
        % @propclass{data} Depends on the kernel setting and problem setup.
        %
        % @default [] @type matrix<double>
        P = [];
    
        % A set of center vectors that is to be used upon calls to
        % evaluateAtCenters.
        %
        % Setting the center values causes some pre-computations to be
        % performed in order for increased performance.
        %
        % @propclass{data}
        %
        % @type matrix<double> @default []
%         Centers;
    end
        
    properties(SetAccess=protected)
        % Flag that determines if the current kernel is a radial basis function, i.e. its
        % evaluation is of the form `\Phi(x,y) = \phi(\noG{x-y})` for some scalar function
        % `\phi`.
        %
        % Set in subclasses according to current kernel.
        %
        % @type logical @default false
        IsRBF = false;
        
        % Flag that determines if the current kernel bases on scalar product evaluations, i.e.
        % are of the form `\Phi(x,y) = \phi(\spG{x}{y})` for some scalar function `\phi`.
        %
        % Set in subclasses according to current kernel.
        %
        % @type logical @default false
        IsScProd = false;
    end
    
    properties(SetAccess=private, GetAccess=protected)
%         fCenters = [];
        fG = 1;
        fP = [];
    end
    
    methods
        function this = BaseKernel
            this = this@KerMorObject;
%             this.registerProps('P', 'G', 'Centers');
            this.registerProps('P', 'G');
        end
        
%         function phix = evaluateAtCenters(this, xi)
%             % Default implementation.
%             %
%             % Returns the kernel matrix for the points `x_i` and the
%             % current Centers
%             phix = this.evaluate(xi, this.fCenters);
%         end
        
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
        
        function copy = clone(this, copy)
            copy.G = this.G;
            copy.P = this.P;
            copy.IsRBF = this.IsRBF;
            copy.IsScProd = this.IsScProd;
        end
    end
    
    %% Getter & Setter
    methods
        function G = get.G(this)
            G = this.fG;
        end
        
        function P = get.P(this)
            P = this.fP;
        end
        
%         function c = get.Centers(this)
%             c = this.fCenters;
%         end
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
        K = evaluate(this, x, y);
        
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

