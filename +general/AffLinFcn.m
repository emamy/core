classdef AffLinFcn < handle
    %AFFLINFCN Summary of this class goes here
    %   Detailed explanation goes here
    % @ingroup general
    
    properties
        Coefficients = cell(1,0);
        Matrices = cell(1,0);
    end
    
    methods
        function res = evaluate(this, t, mu)
            res = 0;
            for idx=1:length(this.Coefficients)
                cfun = this.Coefficients{idx};
                res = res + cfun(t,mu) * this.Matrices{idx};
            end
        end
        
        function addSummand(this, coeff_fcn, mat)
            if ~isa(coeff_fcn,'function_handle')
                error('Parameter "coeff_fcn" must be a function handle');
            end
            if ~isa(mat,'double')
                error('Parameter "mat" must be a double matrix.');
            end
            this.Coefficients{end+1} = coeff_fcn;
            this.Matrices{end+1} = mat;
        end
    end
    
end
