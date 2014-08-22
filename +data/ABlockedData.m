classdef ABlockedData < handle
% ABlockedData: General abstract class that allows computation of and SVD on a large matrix that
% is separated into several blocks.
%
% @author Daniel Wirtz @date 2012-07-11
%
% @change{0,6,dw,2012-10-01} Added an optional Vexclude matrix argument to exclude a given
% space from the data
%
% @new{0,6,dw,2012-07-11} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % The minimum relative value of singular values that triggers selection of the 
        % compared to the largest one.
        %
        % @type double @default 1e-20
        MinRelSingularValueSize = 1e-20;
    end
    
    methods
        function [U, S, V] = getSVD(this, k, Vexclude)
            % Computes an SVD on this blockwise matrix `A = USV^T`.
            %
            % Parameters:
            % k: The number of largest singular values and vectors to compute. @type integer
            % @default all
            % Vexclude: A matrix containing orthonormal columns, whose spanned space is to be
            % excluded from the SVD. @type matrix<double> @default []
            %
            % Return values:
            % U: The left-hand side singular vectors `U` of the decomposition `U\Sigma V^T =
            % A`. @type matrix<double>
            % S: The singular values `\Sigma` of the decomposition `U\Sigma V^T = A`. @type
            % matrix<double>
            % V: The transposed right-hand side singular vectors `V` of the
            % decomposition `U\Sigma V = A`. @type matrix<double>
            [n, m] = size(this);
            if nargin < 3
                Vexclude = [];
                if nargin < 2
                    k = min(n,m);
                end
            end
            % Assume full size if k is empty
            if isempty(k)
                k = min(n,m);
            end
            % Reduce target dimension if remaining dimension is smaller more dimensions are to be excluded
            % than 
            if ~isempty(Vexclude)
                k = min(k, min(n,m)-size(Vexclude,2));
            end
            
            % Fetch case for one block only
            if this.getNumBlocks == 0
                U = [];
                S = [];
            elseif this.getNumBlocks == 1
                Bl = this.getBlock(1);
                % Subtract exclude space if wanted
                if ~isempty(Vexclude)
                    Bl = Bl - Vexclude*(Vexclude'*Bl);
                end
                [U, S] = svd(Bl,'econ');
                U = U(:,1:k);
                S = S(1:k,1:k);
            else            
                % Fully-flavoured singular value decomp using block
                % separation
                opts.issym = 1;
                if KerMor.App.Verbose > 0
                    fprintf('ABlockedData: Computing %d-partial SVD on %dx%d matrix (%d blocks)...\n',...
                        k,n,m,this.getNumBlocks);
                end
                vb = KerMor.App.Verbose > 0 && n*m > 50000 && this.getNumBlocks > 2;
                cnt = 0;
                % If an exclude space is given, choose initial value outside of Vexclude
                if ~isempty(Vexclude)
                    v0 = rand(n,1);
                    opts.v0 = v0 - Vexclude*(Vexclude'*v0);
                end
                doparallel = exist('matlabpool','file') == 2 && matlabpool('size') > 1;
                [U,S] = eigs(@mult,n,k,'la',opts);
                if KerMor.App.Verbose > 2, fprintf('BlockSVD: Finished after %d multiplications.\n',cnt); end
            end
            sel = sqrt(diag(S)/S(1)) >= this.MinRelSingularValueSize & diag(S)>0;
            U = U(:,sel);
            if size(U,2) < k
                warning('KerMor:ABlockedData','Have only %d nonzero singular values instead of %d desired ones.',...
                    size(U,2),k);
            end
            S = sqrt(S(sel,sel));
            if nargout > 2
                V = (U/S)'*this;
            end
            
            function w = mult(v)
                w = 0;
                nb = this.getNumBlocks;
                if vb && ~doparallel
                    pi = ProcessIndicator('ABlockedData: Blockwise multiplication for %d-SVD on %d blocks',...
                        nb,false,k,nb);
                end
                if doparallel
                    if vb
                        fprintf('ABlockedData: Parallel blockwise multiplication for %d-SVD on %d blocks...\n',k,nb);
                    end
                    parfor j = 1:nb
                        B = this.getBlock(j);%#ok
                        % Subtract exclude space if wanted
                        if ~isempty(Vexclude)
                            B = B - Vexclude*(Vexclude'*B);
                        end
                        w = w + B*(B'*v);
                    end
                else               
                    for j = 1:nb
                        B = this.getBlock(j);
                        % Subtract exclude space if wanted
                        if ~isempty(Vexclude)
                            B = B - Vexclude*(Vexclude'*B);
                        end
                        w = w + B*(B'*v);
                        if vb, pi.step; end
                    end
                end
                if vb, pi.stop; end
                cnt=cnt+1;
            end
        end
        
        % Need left-sided matrix multiplication if RHS singular vectors V should be returned.
        function prod = mtimes(matrix, this)          
            if isa(matrix,'data.ABlockedData')
                error('Must override in subclasses for case ABlockedData*matrix.');
            else
                warning('KerMor:ABlockedData',...
                    'If you want the right singular values, you need to override the mtimes method of data.ABlockedData(matrix, this).\nReturning the same instance.');
            end
            prod = this;
        end
        
        function A = toMemoryMatrix(this)
            % Converts this FileMatrix to a full double matrix.
            A = zeros(size(this));
            pos = 1;
            for i=1:this.getNumBlocks
                B = this.getBlock(i);
                s = size(B,2)-1;
                A(:,pos:pos+s) = B;
                pos = pos + s + 1;
            end
        end
    end
    
    methods(Abstract)
        varargout = size(this, dim);
        
        n = getNumBlocks(this);
        
        B = getBlock(this, nr);
    end
end