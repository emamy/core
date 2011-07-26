classdef AffLinCoreFun < dscomponents.ACoreFun & dscomponents.IGlobalLipschitz
%Simple affine-linear core function "f" for a dynamical system.
% 
% Simply wraps an affine-linear function into the ACoreFun interface to
% enable use of simple affine-linear functions as core function. At
% projection, each summand matrix is base changed into the basis given
% by V.
%
% @author Daniel Wirtz @date 15.03.2010
%
% @change{0,5,dw,2011-07-07} Updated this class to use the new general.AffParamMatrix class inside.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        afflinmat;
    end
    
    properties(Dependent, SetAccess=private)
        N;
    end
    
    methods
        function this = AffLinCoreFun
            % Creates a new instance of the AffLinCoreFun.
            this.afflinmat = general.AffParamMatrix;
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            fx = this.afflinmat.compose(t,mu)*x;
        end
        
        function c = getGlobalLipschitz(this, t, mu)
            % Implementation of the interface method from IGlobalLipschitz.
            a = this.afflinmat;
            c = 0;
            for idx=1:length(a.Coefficients)
                cfun = a.Coefficients{idx};
                c = c + abs(cfun(t,mu)) * norm(a.Matrices(:,:,idx));
            end
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            % RHS multiplication of the matrices for correct conversion.
            af = this.afflinmat;
            paf = projected.afflinmat;
            for idx=1:af.N
                paf.Matrices(:,:,idx) = W'*(af.Matrices(:,:,idx)*V);
            end
        end
        
        function addMatrix(this, coeff_fcn, mat)
            this.afflinmat.addMatrix(coeff_fcn, mat);
        end
        
        function N = get.N(this)
            N = this.afflinmat.N;
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = dscomponents.AffLinCoreFun;
            af = general.AffParamMatrix;
            af.Coefficients = this.afflinmat.Coefficients;
            af.Matrices = this.afflinmat.Matrices;
            copy.afflinmat = af;
        end
    end
    
end

