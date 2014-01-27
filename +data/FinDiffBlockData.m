classdef FinDiffBlockData < data.ABlockedData
% FinDiffBlockData: Wrapper for block data that adds finite differences of the block data
%
%
%
% @author Daniel Wirtz @date 2012-09-20
%
% @new{0,6,dw,2012-09-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Access=private)
        orig;
    end
    
    methods
        function this = FinDiffBlockData(original)
            if ~isa(original,'data.ABlockedData')
                error('Only data.ABlockedData instances are allowed');
            end
            this.orig = original;
        end
        
        function prod = mtimes(matrix, this)
            prod = mtimes(matrix,this.orig);
        end
    
        function n = getNumBlocks(this)
            n = this.orig.getNumBlocks;
        end
        
        function B = getBlock(this, nr)
            B = this.orig.getBlock(nr);
            B = [B diff(B,1,2)];
        end
        
        function [n, m] = size(this, dim)
            n = [size(this.orig,1) 2*size(this.orig,2)-this.getNumBlocks];
            if nargin == 2
                if dim > 0 && dim < 3
                    n = n(dim);
                else
                    n = 0;
                end
            elseif nargout == 2
                m = n(2);
                n = n(1);
            end
        end
    end
end