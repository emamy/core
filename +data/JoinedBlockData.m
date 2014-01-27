classdef JoinedBlockData < data.ABlockedData
% JoinedBlockData:
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
        td1;
        td2;
    end
    
    methods
        function this = JoinedBlockData(data1, data2)
            if ~isa(data1,'data.ABlockedData') || ~isa(data2,'data.ABlockedData')
                error('Only joins of two blocked data instances are allowed');
            elseif size(data1,1) ~= size(data2,1)
                error('Can only join blocked data with equal sized first dimensions');
            end
            this.td1 = data1;
            this.td2 = data2;
        end
        
        function prod = mtimes(matrix, this)%#ok
            error('Not yet implemented.');
        end
    
        function n = getNumBlocks(this)
            n = this.td1.getNumBlocks + this.td2.getNumBlocks;
        end
        
        function B = getBlock(this, nr)
            t1len = this.td1.getNumBlocks;
            if nr <= t1len
                B = this.td1.getBlock(nr);
            else
                B = this.td2.getBlock(nr-t1len);
            end
        end
        
        function [n, m] = size(this, dim)
            n = [size(this.td1,1) size(this.td1,2)+size(this.td2,2)];
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