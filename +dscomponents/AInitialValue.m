classdef AInitialValue < dscomponents.IProjectable
% AInitialValue: Abstrat base class for dynamical systems initial values.
%
%
%
% @author Daniel Wirtz @date 2011-07-04
%
% @new{0,5,dw,2011-07-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Abstract)
        % Evaluates the initial value for a given mu.
        %
        % Parameters:
        % mu: The current `\mu` used for simulations, [] otherwise.
        x0 = evaluate(this, mu);
    end
    
end