classdef IRotationInvariant < handle
    % Interface class for rotation invariant kernels
    %
    % All rotation invariant kernels have the form
    % `` \Phi(x) := \phi(||x||), x\in\mathbb{R}^d``
    % for some real-valued scalar function `\phi`.
    %
    % When combinations of Kernels are used, this interface will have to be changed to a property.
    % Up to now, the class CombinationKernel cannot dynamically adopt to the interface for the case
    % that all contained kernels implement this interface.
    %
    % @todo: change property to check for interface implementation,
    % implement in other suitable kernels
    
    methods(Abstract)
        % Allows the evaluation of the function `\phi(x)` for scalar `x`
        % directly.
        Kx = evaluateScalar(x);
    end
    
end

