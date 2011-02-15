classdef AffLinOutputConv < general.AffLinFcn & dscomponents.AOutputConv
    %AFFLINOUTPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function res = project(this, V, W)%#ok
            res = dscomponents.AffLinInputConv;
            res.Coefficients = this.Coefficients;
            n = length(this.Matrices);
            res.Matrices = cell(n,1);
            % RHS multiplication of the matrices for correct conversion.
            for idx=1:n
                res.Matrices{idx} = this.Matrices{idx}*V;
            end
        end
        
        
    end
    
end

