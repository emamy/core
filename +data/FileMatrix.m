classdef FileMatrix < data.FileData & general.ABlockSVD
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
%
% @todo add getBlockPos method for block indices and consider removing the idx vector

    properties(Constant)
        BLOCK_SIZE = 256*1024^2; % 256MB
    end
    
    properties(SetAccess=private)
        % The number of columns for each block.
        %
        % @type integer
        bCols;
        
        % The number of blocks into which this matrix is separated.
        %
        % The number of blocks is computed from the matrix size and the maximum block size in
        % Bytes.
        %
        % @type integer
        %
        % See also: BLOCK_SIZE FileMatrix.FileMatrix
        nBlocks;
        
        % The number of rows
        %
        % @type integer
        n;
        
        % The total number of columns
        %
        % @type integer
        m;
    end
    
    properties(SetAccess=private)
        % An index matrix, storing for each column both the block number and the relative
        % index.
        idx;
        
        % Indicates if a block has explicitly assigned values so far
        created;
        
        % The effective block size
        blocksize;
    end
    
    properties(Access=private, Transient)
        % Cache-related stuff
        cachedBlock = [];
        cachedNr = 0;
        cacheDirty = false;
    end
    
    methods
        function this = FileMatrix(n, m, varargin)
            % Creates a new file matrix.
            %
            % Parameters:
            % n: If a scalar, the row dimension. If a matrix, this must be the only parameter
            % and the file matrix is created using the matrix value and the default BLOCK_SIZE.
            % @type [integer,matrix<double>]
            % m: The number of columns. @type integer
            % varargin: More optional input arguments, see below.
            %
            % Varargin arguments:
            % 1: Either a char array denoting the target root folder where the file-containing
            % data folder should be stored, or a data.FileMatrix instance which is then used to
            % place this FileMatrices files in the same root directory. If the latter is the
            % case, the same block size is also assumed.
            % 2: The maximum block size in Bytes. @type integer @default 256MB
            
            % Matrix case: Create & assign directly
            if nargin < 2
               if ismatrix(n)
                   A = n;
                   [n, m] = size(A);  
               else
                   error('If one argument is passed, it must be a matrix.');
               end
            end
            if length(varargin) < 2
                block_size = min(data.FileMatrix.BLOCK_SIZE,n*m*8);
            elseif isposrealscalar(varargin{2})
                block_size = varargin{2};
            else
                error('Invalid input arguments. Second variable input must be a scalar value determining the block size.');
            end
            if isempty(varargin)
                storage_root = KerMor.App.DataStoreDirectory;
            elseif isa(varargin{1},'data.FileMatrix')
                storage_root = fileparts(varargin{1}.DataDirectory);
                block_size = varargin{1}.blocksize;
            elseif ischar(varargin{1})
                storage_root = varargin{1};
            else
                error('Invalid input arguments. First variable input must be either a FileMatrix or a char array.');
            end
            this = this@data.FileData(fullfile(storage_root,...
                sprintf('matrix_%s',general.IDGenerator.generateID)));
            this.n = n; 
            this.m = m;
            this.blocksize = block_size;
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
            % Overloads the general.ABlockSVD method for efficiency.
            %
            % Transposes this matrix if more rows than columns and multiple blocks are used as
            % it is faster to form the transposed once instead of having a nBlocks^2 cost
            % operation in each call to the custom multiplication function in
            % general.ABlockSVD.getSVD.
            
            if nargin < 2
                k = min(this.n,this.m);
            end
            
            trans = this.nBlocks > 1 && this.n > this.m;
            if trans
                if KerMor.App.Verbose > 0
                    fprintf('FileMatrix.getSVD: Transposing (%dx%d -> %dx%d, %d blocks) for better performance.\n',...
                        this.n,this.m,this.m,this.n,this.nBlocks);
                end
                hlp = this.ctranspose;
            else
                hlp = this;
            end
            if trans || nargout == 3
                [U, S, V] = getSVD@general.ABlockSVD(hlp, k);
            else
                [U, S] = getSVD@general.ABlockSVD(hlp, k);
            end
            % Switch matrices if transposed version was used
            if trans
                tmp = U;
                U = V;
                V = tmp;
            end
        end
        
        function [rowmin, rowmax] = getColBoundingBox(this)
            % Computes the bounding box of the matrix with respect to columns.
            rowmin = Inf*ones(this.n,1);
            rowmax = -rowmin;
            for i=1:this.nBlocks
                A = this.loadBlock(i);
                rowmin = min(rowmin, min(A,[],2));
                rowmax = max(rowmax, max(A,[],2));
            end
        end
        
        function A = toMemoryMatrix(this)
            % Converts this FileMatrix to a full double matrix.
            A = zeros(this.n,this.m);
            for i=1:this.nBlocks
                if this.created(i)
                    pos = (i-1)*this.bCols + 1 : min(i*this.bCols,this.m);
                    A(:,pos) = this.loadBlock(i);
                end
            end
        end
        
        function pos = getBlockPos(this, nr)
            % Returns the column indices of the block "nr" within the full matrix.
            pos = (nr-1)*this.bCols + 1 : min(nr*this.bCols,this.m);
        end
        
        function copy = copyWithNewBlockSize(this, block_size)
            copy = data.FileMatrix(this.n,this.m,fileparts(this.DataDirectory),block_size);
            for k=1:this.nBlocks
                pos = this.getBlockPos(k);
                copy.subsasgn(struct('type',{'()'},'subs',{{':',pos}}),this.loadBlock(k));
            end
        end
        
        %% Overloaded methods
        function s = sum(this, dim)
            if nargin < 2
                dim = 1;
            end
            if dim == 1
                s = zeros(1,this.m);
                for k=1:this.nBlocks
                    s(this.getBlockPos(k)) = sum(this.loadBlock(k),1);
                end
            elseif dim == 2
                s = zeros(this.n,1);
                for k=1:this.nBlocks
                    s = s + sum(this.loadBlock(k),2);
                end
            else
                error('Invalid sum dimension: %d',dim);
            end
        end
        
        function value = power(this, expo)
            if ~isa(this,'data.FileMatrix') && ~isscalar(expo)
                error('FileMatrix power only defined for scalar values.');
            end
            value = data.FileMatrix(this.n,this.m,this);
            for k=1:this.nBlocks
                % Make sure its immediately saved
                value.cacheDirty = true;
                value.saveBlock(k,this.loadBlock(k).^2);
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
                    the_n = this.n;
                    if ~strcmp(s{1},':')
                        the_n = length(s{1});
                    end
                    value = zeros(the_n,length(s{2}));
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
        
        function AB = mtimes(A, B)
            if isa(A,'data.FileMatrix')
                % FileMatrix * FileMatrix case
                if isa(B,'data.FileMatrix')
                    error('Not yet implemented.');
                else
                    if A.m ~= size(B,1)
                        error('Matrix dimensions must agree.');
                    end
                    % FileMatrix * vec case
                    if size(B,2) == 1
                        AB = zeros(A.n,1);
                        for bidx = 1:A.nBlocks
                            Bl = A.loadBlock(bidx);
                            % Only multiply if nonzero
                            if A.created(bidx)
                                AB = AB + Bl*B(A.getBlockPos(bidx),:);  
                            end
                        end
                    % FileMatrix * matrix case
                    else
                        AB = data.FileMatrix(A.n, size(B,2), fileparts(A.DataDirectory),A.blocksize);
                        ABcols = AB.bCols;
                        key = struct('type',{'()'},'subs',{{':'}});
                        for i=1:A.nBlocks
                            if A.created(i)
                                b = A.loadBlock(i);
                                posm = A.getBlockPos(i);
                                % Split up the additions into the block size of AB
                                for j=1:ceil(size(B,2)/ABcols)
                                    key.subs{2} = AB.getBlockPos(j);
                                    % This is reasonably fast as the same block in AB is loaded all
                                    % the time.
                                    AB.subsasgn(key,AB.subsref(key) + b*B(posm,key.subs{2}));
                                end
                            end
                        end
                    end
                end
                % vec * FileMatrix case
            elseif isa(B,'data.FileMatrix')
                AB = data.FileMatrix(size(A,1),B.m, fileparts(B.DataDirectory), B.blocksize);
                key = struct('type',{'()'},'subs',{{':'}});
                for i=1:B.nBlocks
                    if B.created(i)
                        b = B.loadBlock(i);
                        key.subs{2} = B.getBlockPos(i);
                        AB.subsasgn(key,A*b);
                    end
                end
                % Return a matrix if the A argument was a transposed vector
                if size(A,1) == 1
                    AB = AB.toMemoryMatrix;
                end
            end
        end
        
        function trans = ctranspose(this)
            trans = data.FileMatrix(this.m,this.n,...
                fileparts(this.DataDirectory),this.blocksize);
            key = struct('type',{'()'},'subs',{{[],':'}});
            for j=1:this.nBlocks
                B = this.loadBlock(j);
                key.subs{1} = this.getBlockPos(j);
                trans.subsasgn(key,B');
            end
        end
        
        function res = eq(A,B)
            res = false;
            if isa(A,'data.FileMatrix') && all(size(A) == size(B))
                if isa(B,'data.FileMatrix')
                    error('not yet implemented');
                else
                    for i=1:A.nBlocks
                        pos = A.getBlockPos(i);
                        if A.created(i)
                            if ~isequal(A.loadBlock(i),B(:,pos)), return; end
                        else
                            if any(any(B(:,pos) ~= 0)), return; end
                        end
                    end
                    res = true;
                end
            elseif isa(B,'data.FileMatrix')
                res = eq(B, A);
            end
        end
        
        function res = ne(A,B)
            res = ~eq(A,B);
        end
        
        function horzcat(this, varargin)
            error('not yet implemented');
        end
        
        function vertcat(this, varargin)
            error('not yet implemented');
        end
        
        function delete(this)
            if ~this.isSaved
                % Due to caching, no files are written if the entire matrix fits into one
                % block.
                if this.nBlocks > 1
                    % Remove exactly the files for the matrix
                    for k=1:this.nBlocks
                        if this.created(k)
                            delete(fullfile(this.DataDirectory, sprintf('block_%d.mat',k)));
                        end
                    end
                end
            end
            % Superclass delete removes the folder if empty.
            delete@data.FileData(this);
        end
    end
    
    methods(Access=protected)
        function n = getNumBlocks(this)
            n = this.nBlocks;
        end
        
        function n = getColsPerBlock(this)
            n = this.bCols;
        end
        
        function B = getBlock(this, nr)
            B = this.loadBlock(nr);
        end
        
        function [n, m] = getTotalSize(this)
            n = this.n;
            m = this.m;
        end
        
        function this = saveobj(this)
            saveobj@data.FileData(this);
            if this.cacheDirty
                A = this.cachedBlock;%#ok
                save([this.DataDirectory filesep sprintf('block_%d.mat',this.cachedNr)], 'A');
                this.created(this.cachedNr) = true;
                this.cacheDirty = false;
            end
        end
    end
    
    methods(Access=private)
        function saveBlock(this, nr, A)
            % Also save changes to cachedBlock if number matches
            if nr == this.cachedNr
                this.cachedBlock = A;
                this.cacheDirty = true;
            else
                if this.cacheDirty
                    save([this.DataDirectory filesep sprintf('block_%d.mat',nr)], 'A');
                    this.cacheDirty = false;
                end
            end
            this.created(nr) = true;
        end
        
        function A = loadBlock(this, nr)
            if nr == this.cachedNr
                A = this.cachedBlock;
            else
                % Before loading a new block, save the old one if the cache is dirty!
                if this.cacheDirty
                    A = this.cachedBlock;%#ok
                    save([this.DataDirectory filesep sprintf('block_%d.mat',this.cachedNr)], 'A');
                    this.created(this.cachedNr) = true;
                end
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
                this.cachedNr = nr;
                this.cachedBlock = A;
                this.cacheDirty = false;
            end
        end
    end
    
    methods(Static)
        function res = test_FileMatrix
            res = true;
            B = rand(10,100);
            A = data.FileMatrix(10,100,KerMor.App.TempDirectory,1600);
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
            
            % IsEqual test
            res = res && A == B;
            B2 = B; B2(1,:) = -B2(1,:);
            res = res && A ~= B2;
            
            % Transpose test
            At = A';
            res = res && all(all(At.subsref(key) == B'));
            
            % Multiply test
            v = rand(size(A,2),1);
            res = res && isequal(A*v,B*v);
            v = rand(size(A,2),100);
            res = res && A*v == B*v;
            v = rand(1,size(A,1));
            res = res && isequal(v*A,v*B);
            v = rand(100,size(A,1));
            res = res && v*A == v*B;
            
            % SVD test
            p = sqrt(eps);
            [u,s,v] = svd(B,'econ');
            [U,S,V] = A.getSVD;
            V = V.toMemoryMatrix;
            [U5,S5,V5] = A.getSVD(5);
            V5 = V5.toMemoryMatrix;
            res = res && norm(abs(V)-abs(v),'fro') < p && norm(abs(U)-abs(u),'fro') < p &&...
                norm(diag(S)-diag(s)) < p;
            res = res && norm(abs(V5)-abs(v(:,1:5)),'fro') < p ...
                    && norm(abs(U5)-abs(u(:,1:5)),'fro') < p ...
                    && norm(diag(S5)-diag(s(1:5,1:5))) < p;
            
            [ut,st,vt] = svd(B','econ');
            [Ut,St,Vt] = At.getSVD;
            Ut = Ut.toMemoryMatrix;
            [Ut5,St5,Vt5] = At.getSVD(5);
            Ut5 = Ut5.toMemoryMatrix;
            res = res && norm(abs(Vt)-abs(vt),'fro') < p && norm(abs(Ut)-abs(ut),'fro') < p &&...
                norm(diag(St)-diag(st)) < p;
            res = res && norm(abs(Vt5)-abs(vt(:,1:5)),'fro') < p ...
                    && norm(abs(Ut5)-abs(ut(:,1:5)),'fro') < p ...
                    && norm(diag(St5)-diag(st(1:5,1:5))) < p;
                
            % Bounding box test
            [bm, bM] = general.Utils.getBoundingBox(B);
            [am, aM] = A.getColBoundingBox;
            res = res && isequal(bm,am) && isequal(bM,aM);
            
            B = rand(100,10000)+.1;
            A = data.FileMatrix(100,10000,KerMor.App.TempDirectory,round(numel(B)*8/40));
            key.subs = {':', ':'};
            A.subsasgn(key,B);
            
            % Sum test
            bs = sum(B,2);
            as = sum(A,2);
            res = res && isequal(sum(B,1),sum(A,1)) && all(abs((as-bs) ./ bs)) < sqrt(eps);
            
            % Power test
            res = res && A.^2 == B.^2;%#ok
        end
        
        function test_SpeedSVDTransp
            % Test results for 100000x100 matrix:
            % BlockSVD: Computing truncated 50-SVD on 100x1000000 matrix (3 blocks)...
            % Direct time: 2475.53, Transposed time: 756.023, transposed SVD time: 661.171
            %
            % Thus: Auto-transpose for matrices with nBlocks > 1
            A = data.FileMatrix(100,10,KerMor.App.TempDirectory,100*10*8/20);
            A(:,:) = rand(size(A));
            ti = tic;
            [u,s,v] = A.getSVD(5);
            t(1) = toc(ti);
            
            ti = tic;
            B = A';
            t(2) = toc(ti);
            
            ti = tic;
            [v2,s2,u2] = B.getSVD(5);
            t(3) = toc(ti);
            t(2) = t(2) + t(3);
            
            fprintf('Direct time: %g, Transposed time: %g, transposed SVD time: %g\n',t);
        end
    end
end