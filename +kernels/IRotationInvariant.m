classdef IRotationInvariant < handle
    % Interface class for rotation invariant kernels
    %
    % All rotation invariant kernels have the form
    % `` \Phi(x) := \phi(||x||), x\in\mathbb{R}^d``
    % for some real-valued scalar function `\phi`.
    %
    % @todo: change property to check for interface implementation,
    % implement in other suitable kernels
    
    methods(Abstract)
        % Allows the evaluation of the function `\phi(x)` for scalar `x`
        % directly.
        Kx = evaluateScalar(x);
    end
    
end

