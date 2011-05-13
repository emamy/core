classdef IKernelCoeffComp < handle
    % Interface for kernel expansion coefficient computation
    %
    % Any algorithm that can compute onedimensional coefficients
    % given a kernel matrix `K` for a kernel expansion of the type
    % ``f(x) = \sum\limits_{i=1}^N \alpha_i\Phi(x,x_i)``
    % should implement this interface in order to be available as strategy
    % within any approx.BaseCompWiseKernelApprox subclass.
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % @chane{0,3,dw,2011-05-03} Removed offset terms from interface as the `b` offsets for kernel
    % expansions arent used anymore.
    %
    % @new{0,3,dw,2011-03-31} Added this interface.
    
    methods(Abstract)
        % Initialization template method.
        %
        % Call this method before any calls to computeKernelCoefficients
        % and every time the kernel matrix changes.
        %
        % Parameters:
        % K: The Kernel matrix created from the centers `x_i`
        %
        % See also: computeKernelCoefficients
        init(this, K);
        
        % Kernel coefficient computation.
        %
        % Here the concrete class performs the approximation calculation
        % for given function evaluation points `y_i = f(x_i)` at the
        % centers `x_i` also used to compute the kernel matrix passed to 
        % the approx.IKernelCoeffComp#init method.
        %
        % Parameters:
        % yi: The function values `f(x_i)` as row vector.
        %
        % Return values:
        % ai: The coefficients `\alpha_{k,i}` of `f_k(x)`.
        % svidx: The used support vector indices `i` of `x_i`. Optional,
        % leave empty if all are used.
        %
        % See also: init
        [ai, svidx] = computeKernelCoefficients(this, yi);
    end
    
end

