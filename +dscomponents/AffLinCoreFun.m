classdef AffLinCoreFun < dscomponents.ACoreFun ...
        & dscomponents.IGlobalLipschitz & dscomponents.IJacobian
%Simple affine-linear core function "f" for a dynamical system.
% 
% Simply wraps an affine-linear function into the ACoreFun interface to
% enable use of simple affine-linear functions as core function. At
% projection, each summand matrix is base changed into the basis given
% by V.
%
% @author Daniel Wirtz @date 15.03.2010
%
% @new{0,6,dw,2011} Added an optional offset term `b` to the AffLinCoreFun
% to enable affine-linear affine-parametric core functions.
%
% @change{0,5,dw,2011-07-07} Updated this class to use the new general.AffParamMatrix class inside.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=protected)
        AffParamMatrix;
    end
    
    properties
        % The offset term
        b = [];
    end
    
    properties(Dependent, SetAccess=private)
        N;
    end
    
    methods
        function this = AffLinCoreFun
            % Creates a new instance of the AffLinCoreFun.
            this.AffParamMatrix = general.AffParamMatrix;
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            fx = this.AffParamMatrix.compose(t, mu)*x;
            if ~isempty(this.b)
                fx = fx + this.b;
            end
        end
        
        function J = getStateJacobian(this, ~, t, mu)
            J = this.AffParamMatrix.compose(t, mu);
        end
        
        function c = getGlobalLipschitz(this, t, mu)
            % Implementation of the interface method from IGlobalLipschitz.
            a = this.AffParamMatrix;
            c = 0;
            for idx=1:length(a.Coefficients)
                cfun = a.Coefficients{idx};
                c = c + abs(cfun(t,mu)) * norm(a.Matrices(:,:,idx));
            end
        end
        
        function proj = project(this, V, W)
            proj = this.clone;
            proj.AffParamMatrix = this.AffParamMatrix.project(V, W);
            proj.b = W'*this.b;
        end
        
        function addMatrix(this, coeff_fcn, mat)
            this.AffParamMatrix.addMatrix(coeff_fcn, mat);
        end
        
        function N = get.N(this)
            N = this.AffParamMatrix.N;
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            copy = dscomponents.AffLinCoreFun;
            copy.AffParamMatrix = this.AffParamMatrix.clone;
            copy.b = this.b;
        end
    end
    
end

