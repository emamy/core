classdef ManualReduction < spacereduction.BaseSpaceReducer
    % Allows manual selection of the projection matrices `V` and `W`.
    
    properties(Access=private)
        V;
        W;
    end
    
    methods
        
        function this = ManualReduction(V,W)
            % Creates a manual space reducer with given matrices `V` and
            % `W`
            this.V = V;
            if nargin == 1
                if ~isequal(round(V'*V),eye(size(V,2))) || sum(sum(V'*V)-eye(size(V,2))) > 100*eps
                    error('If no W is given, V must be orthogonal!');
                end
                this.W = V;
            else
                this.W = W;
            end
        end
        
        function [V, W] = generateReducedSpace(this, model)%#ok
            % Simply returns the manually selected values.
            V = this.V;
            W = this.W;
        end
        
    end
    
end

