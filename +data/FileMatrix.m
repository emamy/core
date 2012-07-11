classdef FileMatrix < data.FileData
% FileMatrix: File-based matrix which stores sets of rows in separate files.
%
% This class features a caching functionality for the last accessed block, so that for
% subsequent calls to loadBlock with the same number no new hard drive access is necessary.
% This makes the FileMatrix almost as fast as a normal matrix when one block is used, i.e. the
% whole matrix fits into one block of the pre-defined size block_size passed at the
% constructor.
%
% @author Daniel Wirtz @date 2012-07-09
%
% @new{0,6,dw,2012-07-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing   

    properties(Constant)
        BLOCK_SIZE = 256*1024^2; % 256MB
    end
    
    properties
        % The minimum relative value of singular values that triggers selection of the 
        % compared to the largest one.
        MinRelSingularValueSize = 1e-12;
    end

    properties(SetAccess=private)
        n,m;
        
        bCols;
        
        nBlocks;
    end
    
    properties%(Access=private)
        idx;
        
        created;
        
        lastBlock;
        lastNr;
        
        box;
    end
    
    methods
        function this = FileMatrix(n, m, storage_root, block_size)
            % Parameters:
            % block_size: The maximum block size in Bytes. @default 256MB
            
            % Matrix case: Create & assign directly
            if nargin < 2
               if ismatrix(n)
                   A = n;
                   [n, m] = size(A);
               else
                   error('If one argument is passed, it must be a matrix.');
               end
            end
            if nargin < 4
                block_size = min(data.FileMatrix.BLOCK_SIZE,n*m*8);
                if nargin < 3
                   storage_root = KerMor.App.DataStoreDirectory;
                          
                end
            end
            this = this@data.FileData(fullfile(storage_root,...
                sprintf('matrix_%s',general.IDGenerator.generateID)));
            this.n = n; 
            this.m = m;
            this.bCols = max(floor(block_size/(8*n)),1);
            this.nBlocks = ceil(m / this.bCols);
            this.created = false(1,this.nBlocks);
            hlp = reshape(repmat(1:this.nBlocks,this.bCols,1),[],1);
            hlp2 = reshape(repmat(1:this.bCols,this.nBlocks,1)',[],1);
            this.idx = [hlp(1:m) hlp2(1:m)];
            
            % Matrix case: Assign value directly
            if nargin < 2
                this.subsasgn(struct('type',{'()'},'subs',{{':',':'}}),A);
            end
        end
        
        function [U, S, V] = getSVD(this, k)
            psize = min(this.n,this.m);
            if nargin < 2
                k = psize;
            end
            
            opts.issym = 1;
            rowmode = this.m < this.n;
            fprintf('FileMatrix: Computing truncated %d-SVD on %dx%d matrix (%d blocks)...\n',...
                    k,this.n,this.m,this.nBlocks);
            
            fun = @colmode_mult;
            if rowmode
                warning('KerMor:FileMatrix',['Computing SVD on matrix with m < n is very inefficient. '...
                    'Consider arranging the matrix transposed.']);
                fun = @rowmode_mult;    
            end
            [U,S] = eigs(fun,psize,k,'la',opts);
            sel = sqrt(diag(S)/S(1)) >= this.MinRelSingularValueSize;
            U = U(:,sel);
            if size(U,2) < k
                warning('KerMor:FileMatrix','Have only %d nonzero singular values instead of %d desired ones.',...
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
                for j = 1:this.nBlocks
                    Bj = this.loadBlock(j);
                    for i = 1:this.nBlocks
                        Bi = this.loadBlock(i);
                        posi = (i-1)*this.bCols + 1 : min(i*this.bCols,this.m);
                        posj = (j-1)*this.bCols + 1 : min(j*this.bCols,this.m);
                        w(posj,:) = w(posj,:) + Bj'*(Bi*v(posi,:));
                    end
                end
            end
            
            function w = colmode_mult(v)
                w = 0;
                for j = 1:this.nBlocks
                    B = this.loadBlock(j);
                    w = w + B*(B'*v);
                end
            end
        end
        
        function n = numel(~)
            n = 1;
        end
        
        function value = size(this, dim)
            value = [this.n this.m];
            if nargin == 2
                if dim > 0 && dim < 3
                    value = value(dim);
                else
                    value = 1;
                end
            end
        end
        
        function varargout = subsref(this, key)
            % Implements subscripted value retrieval.
            %
            % See also: subsref
            if strcmp(key(1).type,'()')
                s = key.subs;
                % Default case: row & column adressing
                if length(s) == 2
                    % "Select all" mode
                    if strcmp(s{2},':')
                        s{2} = 1:this.m;
                    end
                    pos = this.idx(s{2},:);
                    blocks = unique(pos(:,1));
                    n = this.n;
                    if ~strcmp(s{1},':')
                        n = length(s{1});
                    end
                    value = zeros(n,length(s{2}));
                    for bidx = 1:length(blocks)
                        b = blocks(bidx);
                        B = this.loadBlock(b);
                        value(:,pos(:,1) == b) = B(s{1},pos(pos(:,1)==b,2));
                    end
                % Linear addressing
                elseif length(key.subs) == 1
                    error('Not yet implemented.');
                else
                    error('()-assignment must have one or two subscripts');
                end
                varargout{1} = value;
            else
                [varargout{1:nargout}] = builtin('subsref', this, key);
            end
        end
        
        function this = subsasgn(this, key, value)
            % Implements subscripted assignment.
            %
            % See also: subsasgn
            if strcmp(key(1).type,'()')
                s = key.subs;
                % Default case: row & column adressing
                if length(s) == 2
                    % "Select all" mode
                    if strcmp(s{2},':')
                        s{2} = 1:this.m;
                    end
                    pos = this.idx(s{2},:);
                    blocks = unique(pos(:,1));
                    for bidx = 1:length(blocks)
                        b = blocks(bidx);
                        B = this.loadBlock(b);
                        B(s{1},pos(pos(:,1)==b,2)) = value(:,pos(:,1) == b);
                        this.saveBlock(b,B);
                    end
                % Linear addressing
                elseif length(key.subs) == 1
                    error('Not yet implemented.');
                else
                    error('()-assignment must have one or two subscripts');
                end
            else
                builtin('subsasgn', this, key, value);
            end
        end
        
        function Ax = mtimes(this, x)
            if ~isa(this,'data.FileMatrix') || ~ismatrix(x)
                error('Not yet implemented.');
            elseif this.m ~= size(x,1)
                error('Matrix dimensions must agree.');
            else
                Ax = zeros(this.n,size(x,2));
                for bidx = 1:this.nBlocks
                    B = this.loadBlock(bidx);
                    % Only multiply if nonzero
                    if this.created(bidx)
                        pos = (bidx-1)*this.bCols + 1 : min(bidx*this.bCols,this.m);
                        Ax = Ax + B*x(pos,:);                        
                    end
                end
            end
        end
        
        function trans = ctranspose(this)
            root = fileparts(this.DataDirectory);
            bs = 8*this.bCols*this.n;
            trans = data.FileMatrix(this.m,this.n,root,bs);
            for j=1:this.nBlocks
                B = this.loadBlock(j);
                pos = (j-1)*this.bCols + 1 : min(j*this.bCols,this.m);
                trans.subsasgn(struct('type',{'()'},'subs',{{pos,':'}}),B');
            end
        end
        
        function delete(this)
            if ~this.isSaved
                % Remove exactly the files for the matrix
                for k=1:this.nBlocks
                    if this.created(k)
                        delete(fullfile(this.DataDirectory, sprintf('block_%d.mat',k)));
                    end
                end
            end
            % Superclass delete removes the folder if empty.
            delete@data.FileData(this);
        end
        
        function [rowmin, rowmax] = getColBoundingBox(this)
            rowmin = this.box(:,1);
            rowmax = this.box(:,2);
        end
    end
    
    methods(Access=private)
        function saveBlock(this, nr, A)
            % update bounding box
            if isempty(this.box)
                this.box(:,1) = min(A,[],2);
                this.box(:,2) = max(A,[],2);
            else
                this.box(:,1) = min([this.box(:,1) min(A,[],2)],[],2);
                this.box(:,2) = min([this.box(:,2) max(A,[],2)],[],2);
            end
            save([this.DataDirectory filesep sprintf('block_%d.mat',nr)], 'A');
            % Also save changes to lastBlock if number matches
            if isequal(nr, this.lastNr)
                this.lastBlock = A;
            end
            this.created(nr) = true;
        end
        
        function A = loadBlock(this, nr)
            if isequal(nr, this.lastNr)
                A = this.lastBlock;
            else
                if this.created(nr)
                    s = load([this.DataDirectory filesep sprintf('block_%d.mat',nr)]);
                    A = s.A;
                else
                    A = zeros(this.n, this.bCols);
                    % Shorten the last block to effectively used size
                    if nr == this.nBlocks
                        A = A(:,1:(this.m-(this.nBlocks-1)*this.bCols));
                    end
                end
                this.lastNr = nr;
                this.lastBlock = A;
            end
        end
    end
    
    methods(Static)
        function res = test_FileMatrix
            res = true;
            B = rand(10,100);
            A = data.FileMatrix(10,100,KerMor.App.DataStoreDirectory,1600);
            key = struct('type',{'()'},'subs',[]);
            % col-wise setting (fast)
            for k=1:100
                key.subs = {':', k};
                A.subsasgn(key,B(:,k));
            end
            key.subs = {':', ':'};
            res = res && all(all(A.subsref(key) == B));
            % Row-wise setting
            for k=1:10
                key.subs = {k,':'};
                A.subsasgn(key,B(k,:));
            end
            key.subs = {':', ':'};
            res = res && all(all(A.subsref(key) == B));
            
            % Direct constructor test
            A = data.FileMatrix(B);
            res = res && all(all(A.subsref(key) == B));
            
            % Transpose test
            At = A';
            res = res && all(all(At.subsref(key) == B'));
            
            % Multiply test
            v = rand(size(A,2),1);
            res = res && isequal(A*v,B*v);
            v = rand(size(A,2),100);
            res = res && isequal(A*v,B*v);
            
            % SVD test
            p = sqrt(eps);
            [u,s,v] = svd(B,'econ');
            [U,S,V] = A.getSVD;
            [U5,S5,V5] = A.getSVD(5);
            res = res && norm(abs(V)-abs(v),'fro') < p && norm(abs(U)-abs(u),'fro') < p &&...
                norm(diag(S)-diag(s)) < p;
            res = res && norm(abs(V5)-abs(v(:,1:5)),'fro') < p ...
                    && norm(abs(U5)-abs(u(:,1:5)),'fro') < p ...
                    && norm(diag(S5)-diag(s(1:5,1:5))) < p;
            
            [ut,st,vt] = svd(B','econ');
            % Arguments exchanged here as first one is always the smaller matrix
            [Ut,St,Vt] = At.getSVD;
            [Ut5,St5,Vt5] = At.getSVD(5);
            res = res && norm(abs(Vt)-abs(vt),'fro') < p && norm(abs(Ut)-abs(ut),'fro') < p &&...
                norm(diag(St)-diag(st)) < p;
            res = res && norm(abs(Vt5)-abs(vt(:,1:5)),'fro') < p ...
                    && norm(abs(Ut5)-abs(ut(:,1:5)),'fro') < p ...
                    && norm(diag(St5)-diag(st(1:5,1:5))) < p;
                
            % Bounding box test
            [bm, bM] = general.Utils.getBoundingBox(B);
            [am, aM] = A.getColBoundingBox;
            res = res && isequal(bm,am) && isequal(bM,aM);
        end
    end
end