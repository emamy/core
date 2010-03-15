classdef AffLinCoreFcn < dscomponents.ICoreFun
    %Simple affine-linear core function "f" for a dynamical system.
    % 
    % Simply wraps an affine-linear function into the ICoreFun interface to
    % enable use of simple affine-linear functions as core function. At
    % projection, each summand matrix is base changed into the basis given
    % by V.
    %
    % @Daniel Wirtz, 15.03.2010
    
    properties(Access=private)
        afflinfcn = general.AffLinFcn;
    end
    
    methods
        function fx = evaluate(this, x, t, mu)
            fx = this.afflinfcn.evaluate(t,mu)*x;
        end
        
        function res = project(this, V)
            res = AffLinCoreFcn;
            af = general.AffLinFcn;
            af.Coefficients = this.afflinfcn.Coefficients;
            af.Matrices = cell(n,1);
            % RHS multiplication of the matrices for correct conversion.
            for idx=1:n
                ad.Matrices{idx} = V'*this.afflinfcn.Matrices{idx}*V;
            end
            res.afflinfcn = af;
        end
        
        function addSummand(this, coeff_fcn, mat)
            this.afflinfcn.addSummand(coeff_fcn, mat);
        end
    end
    
end

