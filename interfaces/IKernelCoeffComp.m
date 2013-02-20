classdef IKernelCoeffComp < handle
    % Interface for kernel expansion coefficient computation
    %
    % Any algorithm that can compute onedimensional coefficients
    % given a kernel matrix `K` for a kernel expansion of the type
    % ``f(x) = \sum\limits_{i=1}^N \alpha_i\Phi(x,x_i)``
    % should implement this interface in order to be available as strategy
    % within any approx.KernelApprox subclass.
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % @change{0,5,dw,2011-09-12}
    % - Added initial values that can be passed to the CoeffComp algorithms.
    % - New Property
    % IKernelCoeffComp.MultiTargetComputation that
    % indicates of the algorithm at hand can handle a column vector matrix
    % instead of a single vector.
    %
    % @new{0,5,dw,2011-07-07} Moved this class to the approx.algorithms package.
    %
    % @change{0,3,dw,2011-05-03} Removed offset terms from interface as the `b` offsets for kernel
    % expansions arent used anymore.
    %
    % @new{0,3,dw,2011-03-31} Added this interface.
    
    properties(SetAccess=protected)
        % A flag that indicates to users if the coefficient computation
        % method is capable of using a matrix of column fxi vectors or only
        % single vectors.
        %
        % @type logical @default false
        MultiTargetComputation = false;
    end
    
    methods(Abstract)
        % Initialization template method.
        %
        % Call this method before any calls to computeKernelCoefficients
        % and every time the kernel matrix changes.
        %
        % Parameters:
        % K: The Kernel matrix for the centers `x_i` @type data.FileMatrix
        %
        % See also: computeKernelCoefficients
        init(this, K);
        
        % Kernel coefficient computation.
        %
        % Here the concrete class performs the approximation calculation
        % for given function evaluation points `y_i = f(x_i)` at the
        % centers `x_i` also used to compute the kernel matrix passed to 
        % the IKernelCoeffComp#init method.
        %
        % Parameters:
        % yi: The function values `f(x_i)` as column vector. If
        % MultiTargetComputation is true, this can be a matrix of column
        % vectors.
        % initialai: The values to use as initial coefficients. It is up to
        % the implementing classes to use those if passed; however, a call
        % with empty argument must be possible, too.
        %
        % Return values:
        % ci: The coefficients `c_{k,i}` of `f_k(x)` as ROW vector. @type rowvec
        % svidx: The used support vector indices `i` of `x_i`. Always
        % required. (Set to 1:n if no sparsity is given by the coeffcomp
        % method) @type integer
        %
        % See also: init
        [ci, svidx] = computeKernelCoefficients(this, yi, initialai);
    end
    
    methods(Abstract, Static)
        % Returns a default set of configurations for this IKernelCoeffComp
        c = getDefaultConfig;
    end
    
end

