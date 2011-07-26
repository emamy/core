classdef AffLinOutputConv < general.AffParamMatrix & dscomponents.AOutputConv
    %AFFLINOUTPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function C = evaluate(this, t, mu)
            C = this.compose(t,mu);
        end
        
        function res = project(this, V, W)%#ok
            projected = this.clone;
            % RHS multiplication of the matrices for correct conversion.
            for idx=1:this.N
                projected.Matrices(:,:,idx) = this.Matrices(:,:,idx)*V;
            end
        end
        
        function copy = clone(this)
            copy = dscomponents.AffLinOutputConv;
            copy = clone@general.AffParamMatrix(this, copy);
        end
    end
    
end

