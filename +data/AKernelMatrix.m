classdef AKernelMatrix < KerMorObject
% AKernelMatrix: 
%
% Abstract base class for KernelMatrices, either memory or file-based.
%
% @author Daniel Wirtz @date 2011-09-09
%
% @new{0,5,dw,2011-09-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Abstract)
        % Subscript access overload
        value = subsref(this, key);
        
        % Subscript assignment overload
        this = subsasgn(this, key, value)
        
        % An multiplication overload
        prod = mtimes(A, B)
        
        % Augments the kernel matrix by the given vector, whos length must
        % be one bigger than the current matrices size.
        augment(this, vec)
    end
    
end