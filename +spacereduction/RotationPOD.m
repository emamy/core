classdef RotationPOD < spacereduction.PODReducer
% RotationPOD: Implements the normal PODReducer but rotates the resulting matrix.
%
% This space-reduction class has been added in order to be able to obtain slightly worse subspace
% computations as would be given by a full POD.
%
% @author Daniel Wirtz @date 2011-05-23
%
% @new{0,4,dw,2011-05-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % `\in [0, 2\pi]`
        %
        % @default .1
        Degree = pi/.003;
        
        % Number of dimensions to rotate by degree
        Dims = 5;
        
        R;
    end
    
    methods
        function [V,W] = generateReducedSpace(this, model)
            V = generateReducedSpace@spacereduction.PODReducer(this, model);
            
            R = eye(size(V,1)); %#ok<*PROP>
            for idx = 1:this.Dims
                Q = eye(size(V,1));
                idx1 = randi(this.Dims);
                idx2 = randi(this.Dims);
                Q(idx1,idx1) = cos(this.Degree);
                Q(idx1,idx2) = -sin(this.Degree);
                Q(idx2,idx1) = sin(this.Degree);
                Q(idx2,idx2) = cos(this.Degree);
                R = R*Q;
            end
            V = R*V;
            W = V;
            %this.R = R;
        end
    end
    
end