classdef AffLinOutputConv < general.AffParamMatrix & dscomponents.AOutputConv
    %AFFLINOUTPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function C = evaluate(this, t, mu)
            C = this.compose(t,mu);
        end
        
        function projected = project(this, V, W)
            projected = project@general.AProjectable(this, V, W, this.clone);
            % RHS multiplication of the matrices for correct conversion.
            projected.Matrices = zeros(size(V,2),this.N);
            for idx=1:this.N
                projected.Matrices(:,idx) = reshape(this.getMatrix(idx)*V,[],1);
            end
        end
        
        function copy = clone(this)
            copy = clone@general.AffParamMatrix(this, ...
                dscomponents.AffLinOutputConv);
        end
    end
    
end

