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
        function [U, S, V] = getSVD(this, k, Vexclude, targetdims)
            % Computes an SVD on this blockwise matrix `A = USV^T`.
            %
            % Parameters:
            % k: The number of largest singular values and vectors to compute. @type integer
            % @default all
            % Vexclude: A matrix containing orthonormal columns, whose
            % spanned space is to be excluded from the SVD. If targetdims
            % are given, the first dimension must match the targeted
            % dimension's size. @type matrix<double> @default []
            % targetdims: The dimensions on which to perform the SVD.
            % Selects all by default. @type colvec<integer> @default ':'
            %
            % Return values:
            % U: The left-hand side singular vectors `U` of the decomposition `U\Sigma V^T =
            % A`. @type matrix<double>
            % S: The singular values `\Sigma` of the decomposition `U\Sigma V^T = A`. @type
            % matrix<double>
            % V: The transposed right-hand side singular vectors `V` of the
            % decomposition `U\Sigma V = A`. @type matrix<double>
            if nargin < 4
                targetdims = ':';
                if nargin < 3
                    Vexclude = [];
                end
            end
            % Get the correct dimensions if some are to be excluded
            [n, m] = size(this);
            if ~ischar(targetdims)
                % Use the number of target dims for n if not all are used
                n = length(targetdims);
            end
            % Assume full size if k is not set or empty
            if nargin < 2 || isempty(k)
                k = min(n,m);
            end
            % Reduce target dimension if remaining dimension is smaller
            if ~isempty(Vexclude)
                k = min(k, min(n,m)-size(Vexclude,2));
            else
                k = min(k, min(n,m));
            end
            
            % Fetch case for one block only
            if this.getNumBlocks == 0
                U = [];
                S = [];
            elseif this.getNumBlocks == 1
                Bl = this.getBlock(1);
                Bl = Bl(targetdims,:);
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
                vb = KerMor.App.Verbose > 1 && n*m > 50000 && this.getNumBlocks > 2;
                cnt = 0;
                % If an exclude space is given, choose initial value outside of Vexclude
                if ~isempty(Vexclude)
                    v0 = rand(n,1);
                    opts.v0 = v0 - Vexclude*(Vexclude'*v0);
                end
                doparallel = exist('matlabpool','file') == 2 && matlabpool('size') > 1;
                [U,S] = eigs(@mult,n,k,'la',opts);
                if KerMor.App.Verbose > 2, fprintf('BlockSVD: Finished after %d multiplications.\n',cnt); end
                S = sqrt(S);
            end
            sel = sqrt(diag(S)/S(1)) >= this.MinRelSingularValueSize & diag(S)>0;
            U = U(:,sel);
            if size(U,2) < k
                warning('KerMor:ABlockedData','Have only %d nonzero singular values instead of %d desired ones.',...
                    size(U,2),k);
            end
            S = S(sel,sel);
            if nargout > 2
                % Treat this special case here instead of cloning everything
                if ~ischar(targetdims)
                    error('Fixme: Not correctly implemented yet.');
                    [nfull,~] = size(this);
                    mat = eye(nfull);
                    mat(targetdims,targetdims) = (U/S)';
                    V = mat*this;
                    V = V(targetdims,:);
                else
                    V = (U/S)'*this;
                end
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
                        B = B(targetdims,:);
                        % Subtract exclude space if wanted
                        if ~isempty(Vexclude)
                            B = B - Vexclude*(Vexclude'*B);
                        end
                        w = w + B*(B'*v);
                    end
                else               
                    for j = 1:nb
                        B = this.getBlock(j);
                        B = B(targetdims,:);
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
                    'If you want the right singular values, you need to override the mtimes method of data.ABlockedData(matrix, this).\nReturning the same instance for V.');
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
    
    methods(Static)
        function res = test_BlockSVD_vs_SVD
            blocks = 5;
            ns = 6;
            an = [500 1000 700];
            am = [500 700 1000];
            res = true;
            for v = 1:3
                n = an(v);
                m = am(v);
                for nb=1:blocks
                    B = rand(n,m);
                    A = data.FileMatrix(n,m,'BlockSize',data.FileMatrix.blockSizeOf(B)/nb);
                    A.subsasgn(struct('type',{'()'},'subs',{{':', ':'}}),B);
                    [ub,us] = A.getSVD(ns);
                    [u,s] = svd(B,'econ');
                    chk = norm(diag(us)-diag(s(1:ns,1:ns)));
                    res = res && chk < 1e-8;
                    % mess with directions - so correct
                    ubs = sign(ub(1,:));
                    us = sign(u(1,1:ns));
                    chk2 = norm(bsxfun(@times,ub,ubs)-bsxfun(@times,u(:,1:ns),us));
                    res = res & chk2 < 1e-8;
                    fprintf('Comparison Block-SVD vs SVD (%d singular values): Singular value vector norm diff: %g, U norm diff: %g\n',ns,chk,chk2);
                end
            end
        end
    end
end