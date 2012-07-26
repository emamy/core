classdef ABlockSVD < handle
% ABlockSVD: General abstract class that allows computation of and SVD on an object that is
% separated into several blocks.
%
% This implementation favors storage of the vectors column-wise, i.e. SVDs on matrices
% `A\in\R^{n\times m}` with `m\geq n` are more efficient as each block has to be loaded only
% once.
%
% @author Daniel Wirtz @date 2012-07-11
%
% @new{0,6,dw,2012-07-11} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The minimum relative value of singular values that triggers selection of the 
        % compared to the largest one.
        %
        % @type double @default 1e-12
        MinRelSingularValueSize = 1e-12;
    end
    
    methods
        function [U, S, V] = getSVD(this, k)
            % Computes an SVD on this block matrix `A`.
            %
            % Parameters:
            % k: The number of largest singular values and vectors to compute. @type integer
            % @default all
            %
            % Return values:
            % U: The left-hand side singular vectors `U` of the decomposition `U\Sigma V^T =
            % A`. @type matrix<double>
            % S: The singular values `\Sigma` of the decomposition `U\Sigma V^T = A`. @type
            % matrix<double>
            % V: The right-hand side singular vectors `V` of the decomposition `U\Sigma V^T =
            % A`. @type matrix<double>
            [n, m] = this.getTotalSize;
            psize = min(n,m);
            if nargin < 2
                k = psize;
            end
            
            opts.issym = 1;
            rowmode = m < n;
            if KerMor.App.Verbose > 0
                fprintf('BlockSVD: Computing truncated %d-SVD on %dx%d matrix (%d blocks)...\n',...
                    k,n,m,this.getNumBlocks);
            end
            
            fun = @colmode_mult;
            if rowmode
                fun = @rowmode_mult;    
                if this.getNumBlocks > 1
                    warning('KerMor:BlockSVD',['Computing SVD on matrix with m < n is very inefficient. '...
                        'Consider arranging the matrix transposed.']);
                end
            end
            [U,S] = eigs(fun,psize,k,'la',opts);
            sel = sqrt(diag(S)/S(1)) >= this.MinRelSingularValueSize;
            U = U(:,sel);
            if size(U,2) < k
                warning('KerMor:BlockSVD','Have only %d nonzero singular values instead of %d desired ones.',...
                    size(U,2),k);
            end
            S = sqrt(S(sel,sel));
            if rowmode
                hlp = this*(U/S);
                V = U;
                U = hlp;
            elseif nargout > 2
                V = this'*(U/S);
            end
            
            function w = rowmode_mult(v)
                w = zeros(size(v));
                nb = this.getNumBlocks;
                nc = this.getColsPerBlock;
                for j = 1:nb
                    Bj = this.getBlock(j);
                    for i = 1:nb
                        Bi = this.getBlock(i);
                        posi = (i-1)*nc + 1 : min(i*nc,m);
                        posj = (j-1)*nc + 1 : min(j*nc,m);
                        w(posj,:) = w(posj,:) + Bj'*(Bi*v(posi,:));
                    end
                end
            end
            
            function w = colmode_mult(v)
                w = 0;
                nb = this.getNumBlocks;
                pi = tools.ProcessIndicator('ABlockSVD: Blockwise multiplication for SVD on %d blocks',nb,...
                    false,nb);
                for j = 1:nb
                    B = this.getBlock(j);
                    w = w + B*(B'*v);
                    pi.step;
                end
                pi.stop;
            end
        end
    end
    
    methods(Abstract, Access=protected)
        n = getNumBlocks(this);
        
        n = getColsPerBlock(this);
        
        B = getBlock(this, nr);
        
        [n, m] = getTotalSize(this);
    end
    
end