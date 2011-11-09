classdef IKernelMatrix < KerMorObject & ICloneable
% IKernelMatrix: Abstract base interface for KernelMatrices, either memory
% or file-based.
%
% @author Daniel Wirtz @date 2011-09-09
%
% @new{0,5,dw,2011-09-09} New template method plus to force implementation of addition.
%
% @new{0,5,dw,2011-09-09} 
% - New template method mldivide to force implementation of matrix left division
% - Implementing ICloneable now
%
% @new{0,5,dw,2011-09-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetAccess=protected)
        % The dimension of the kernel matrix.
        %
        % @type integer
        Dim = [];
    end
    
    methods(Abstract)
        % Subscript access overload
        value = subsref(this, key);
        
        % Subscript assignment overload
        this = subsasgn(this, key, value)
        
        % An multiplication overload
        prod = mtimes(A, B)
        
        % An addition overload
        prod = plus(A, B)
        
        % Augments the kernel matrix by the given vector, whos length must
        % be one bigger than the current matrices size.
        augment(this, vec)
        
        % Left division `K\backslash f`
        a = mldivide(this, vec)
    end
    
end