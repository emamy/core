classdef AffLinOutputConv < general.AffParamMatrix & dscomponents.AOutputConv
    %AFFLINOUTPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function C = evaluate(this, t, mu)
            C = this.compose(t,mu);
        end
        
        function projected = project(this, V, W)%#ok
            % Dont store V,W due to hard drive space saving (not really needed here)
            %projected = project@general.AProjectable(this, V, W, this.clone);
            
            projected = this.clone;
            % RHS multiplication of the matrices for correct conversion.
            for idx=1:this.N
                projected.Matrices(:,:,idx) = this.Matrices(:,:,idx)*V;
            end
        end
        
        function copy = clone(this)
            copy = clone@general.AffParamMatrix(this, ...
                dscomponents.AffLinOutputConv);
        end
    end
    
end

