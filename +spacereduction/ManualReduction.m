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
                if ~isequal(round(V'*V),eye(size(V,2))) || sum(sum(V'*V-eye(size(V,2)))) > 100*eps
                    error('If no W is given, V must be orthogonal!');
                end
                this.W = V;
            else
                if ~isequal(round(W'*V),eye(size(W,2))) || sum(sum(W'*V-eye(size(W,2)))) > sqrt(eps)
                    error('V and W must be bi-orthogonal!');
                end
                this.W = W;
            end
        end
        
    end
    
    methods(Access=protected)
        
        function [V,W] = generateReducedSpaceImpl(this, model, subset)%#ok
            % Simply returns the manually selected values.
            V = this.V(subset,:);
            W = [];
            if ~isempty(this.W)
                W = this.W(subset,:);
            end
        end
        
    end
    
end

