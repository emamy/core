classdef PointerInitialValue < dscomponents.AInitialValue
% PointerInitialValue: Allows initial values using function pointers for actual evaluation.
%
% @author Daniel Wirtz @date 2011-11-15
%
% @new{0,6,dw,2011-11-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        target;
    end
    
    methods
        function this = PointerInitialValue(ptr)
            % Creates a new pointer initial value.
            %
            % Parameters:
            % ptr: A function handle taking mu as argument @type function_handle
            this = this@dscomponents.AInitialValue;
            if isa(ptr,'function_handle')
                this.target = ptr;
            else
                error('Argument ptr must be a function handle.');
            end
        end
        
        function x0 = evaluate(this, mu)
            % Evaluates the initial value for a given mu.
            %
            % Parameters:
            % mu: The current `\mu` used for simulations, [] otherwise.
            x0 = this.target(mu);
        end
        
        function proj = project(this, V, W)%#ok
            % Projects the initial value function into the reduced space.
            % Creates a new PointerInitialValue and computes `W^tx_0(\mu)`
            
            newfun = @(mu)W' * this.target(mu);
            proj = dscomponents.PointerInitialValue(newfun);
            % Dont store V,W due to hard drive space saving (not really needed here, W will be
            % stored in function handle anyways)
            %proj = project@general.AProjectable(this, V, W, proj);
        end
        
        function copy = clone(this)
            copy = dscomponents.PointerInitialValue(this.target);
        end
    end
    
end