classdef AxBlockData < data.ABlockedData
% AxBlockData: Wrapper for block data that computes A(x) of the block data
%
% No time dependency allowed yet.
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
        fun;
        
        % The scaled model times for each trajectory
        times;
    end
    
    methods
        function this = AxBlockData(original, fun, times)
            if ~isa(original,'data.ATrajectoryData')
                error('Only data.ATrajectoryData instances can be wrapped');
            end
            if ~isa(fun, 'dscomponents.ACoreFun')
                error('Fun argument must be a dscomponents.ACoreFun');
            end
            this.orig = original;
            this.fun = fun;
            this.times = times;
        end
        
        function prod = mtimes(matrix, this)
            prod = mtimes(matrix,this.orig);
        end
    
        function n = getNumBlocks(this)
            n = this.orig.getNumBlocks;
        end
        
        function B = getBlock(this, nr)
            [B, mu] = this.orig.getTrajectoryNr(nr);
            if ~this.fun.TimeDependent
                B = this.fun.evaluate(B,[],mu);
            else
                B = this.fun.evaluate(B,this.times,repmat(mu,1,length(this.times)));
            end
        end
        
        function [n, m] = size(this, dim)
            n = size(this.orig);
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