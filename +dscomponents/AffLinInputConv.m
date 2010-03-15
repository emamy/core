classdef AffLinInputConv < general.AffLinFcn & dscomponents.IInputConv
    %AFFLININPUTCONV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function res = project(this, V)
            res = dscomponents.AffLinInputConv;
            res.Coefficients = this.Coefficients;
            n = length(this.Matrices);
            res.Matrices = cell(n,1);
            
            % LHS multiplication of the matrices for correct conversion.
            for idx=1:n
                res.Matrices{idx} = V'*this.Matrices{idx};
            end
        end
        
        
    end
    
end

