classdef AffLinFcn < KerMorObject
    %AFFLINFCN Summary of this class goes here
    %   Detailed explanation goes here
    % @ingroup general
    %
    % @change{0,3,sa,2011-05-06} Implemented Setters for the class
    % properties
    
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
        
        function set.Coefficients(this, value)
            if ~iscell(value)
                error('Value must be a cell array');
            end
            this.Coefficients = value;
        end
        
        function set.Matrices(this, value)
            if ~iscell(value)
                error('Value must be a cell array');
            end
            this.Matrices = value;
        end
        
        function addSummand(this, mat)
           
            if ~isa(mat,'double')
                error('Parameter "mat" must be a double matrix.');
            end
            
            this.Matrices{end+1} = mat;
        end
    end
    
end

