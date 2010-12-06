classdef AffLinCoreFun < dscomponents.ACoreFun & dscomponents.IGlobalLipschitz
    %Simple affine-linear core function "f" for a dynamical system.
    % 
    % Simply wraps an affine-linear function into the ACoreFun interface to
    % enable use of simple affine-linear functions as core function. At
    % projection, each summand matrix is base changed into the basis given
    % by V.
    %
    % @author Daniel Wirtz @date 15.03.2010
    
    properties(Access=private)
        afflinfcn;
    end
    
    methods
        function this = AffLinCoreFun
            % Creates a new instance of the AffLinCoreFun.
            this.afflinfcn = general.AffLinFcn;
        end
        
        function updateSimConstants(this)%#ok
            % Nothing to do here.
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            fx = this.afflinfcn.evaluate(t,mu)*x;
        end
        
        function c = getGlobalLipschitz(this, t, mu)
            % Implementation of the interface method from IGlobalLipschitz.
            a = this.afflinfcn;
            c = 0;
            for idx=1:length(a.Coefficients)
                cfun = a.Coefficients{idx};
                c = c + abs(cfun(t,mu)) * norm(a.Matrices{idx});
            end
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            % RHS multiplication of the matrices for correct conversion.
            for idx=1:length(this.afflinfcn.Matrices)
                projected.afflinfcn.Matrices{idx} = W'*this.afflinfcn.Matrices{idx}*V;
            end
        end
        
        function addSummand(this, coeff_fcn, mat)
            this.afflinfcn.addSummand(coeff_fcn, mat);
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = dscomponents.AffLinCoreFun;
            af = general.AffLinFcn;
            af.Coefficients = this.afflinfcn.Coefficients;
            af.Matrices = this.afflinfcn.Matrices;
            copy.afflinfcn = af;
        end
    end
    
end

