classdef FileMatrix < data.FileData & data.ABlockedData
    % FileMatrix: File-based matrix which stores sets of columns in separate files.
    %
    % This class features a caching functionality for the last accessed block, so that for
    % subsequent calls to loadBlock with the same number no new hard drive access is necessary.
    % This makes the FileMatrix almost as fast as a normal matrix when one block is used, i.e. the
    % whole matrix fits into one block of the pre-defined size block_size passed at the
    % constructor.
    %
    % @author Daniel Wirtz @date 2012-07-09
    %
    % @change{0,6,dw,2012-11-06} Changed the constructor to use an inputParser
    %
    % @new{0,6,dw,2012-07-09} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    %
    % @todo add getBlockPos method for block indices and consider removing the idx vector
    
    properties(Constant)
        % The default block size to use for new FileMatrix instances.
        %
        % This is used as default value in KerMor.App.BlockSize, if not set
        % differently (machine dependent)
        %
        % The block size unit is Megabyte.
        %
        % @type integer @default 512MB
        BLOCK_SIZE = 512; %[MB]
    end
    
    properties
        % Set this flag to true, if only so much space should be required as actually needed
        % when transposing the matrix (i.e. twice the size of the current matrix).
        %
        % If false (default), running ctranspose will create small blocks fitted to the size of
        % the meshed sizes temporarily. This needs 3 times the space instead of 2 but avoids
        % reloading the whole block even if no data is changed where it already exists.
        InPlaceTranspose = false;
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
        
        MinValue = Inf;
        
        MaxValue = -Inf;
    end
    
    properties(Access=private)
        % An index matrix, storing for each column both the block number and the relative
        % index.
        idx;
        
        % Indicates if a block has explicitly assigned values so far
        created;
        
        % The effective block size
        blocksize;
    end
    
    properties(Access=private)
        % Cache-related stuff
        cachedBlock = [];
        cachedNr = 0;
        cacheDirty = false;
    end
    
    methods
        function this = FileMatrix(var, varargin)
            % Creates a new file matrix.
            % Possible constructors:
            % - FileMatrix(A): Creates a new file matrix
            %
            % Parameters:
            % var: If a scalar, the row dimension. If a matrix, the file matrix is initialized
            % using the matrix value. @type matrix<double>
            % varargin: More optional input arguments, see below.
            % m: If the var argument was a row dimension, the column dimension `m` is a
            % required argument. @type integer
            % Dir: A char array denoting the target root folder where the file-containing
            % data folder should be stored.
            % BlockSize: The maximum block size in MegaBytes. @type integer
            
            if isempty(var)
                error('Cannot create a FileMatrix with empty first argument. Must either be a row number or a matrix.');
            end
            
            ip = inputParser;
            matrixin = false;
            if numel(var) > 1
                A = var;
                [n, m] = size(A); %#ok<*PROP>
                matrixin = true;
%             elseif isa(var,'data.FileMatrix')
%                 n = var.n;
%                 m = var.m;
%                 varargin = {'Dir', fileparts(var.DataDirectory),...
%                     'BlockSize', var.blocksize};
            else
                n = var;
                ip.addRequired('m');
            end
            
            ip.addParameter('Dir',KerMor.App.TempDirectory,...
                @(v)isempty(v) || (ischar(v) && exist(v,'dir') == 7));
            ip.addParameter('BlockSize',KerMor.App.BlockSize,@(v)isposrealscalar(v));
            ip.parse(varargin{:});
            % Latest at here we have an m value
            if isfield(ip.Results,'m')
                m = ip.Results.m;
            end
            % Error checks
            if ~(round(n) == n && round(m) == m)
                error('Size arguments n, m must be integer values.');
            end
            BYTEPERDOUBLE = 8;
            bsize = min(ip.Results.BlockSize*1024^2,n*m*BYTEPERDOUBLE);
            bcols = max(floor(bsize/(BYTEPERDOUBLE*n)),1);
            nb = ceil(m / bcols);
            
            % Determine the DataDirectory (store in each case, even if there is only one block.
            % This is due to the fact that matrix operations can result in a filematrix with
            % more than one block, thus needing to know where the previous matrix should have
            % been located.
            ddir = ip.Results.Dir;
            if isempty(ddir)
                ddir = KerMor.App.TempDirectory;
            end
            this = this@data.FileData(fullfile(ddir,...
                sprintf('matrix_%s',IDGenerator.generateID)));
            
            this.n = n;
            this.m = m;
            this.blocksize = bsize/1024^2;
            this.bCols = bcols;
            this.nBlocks = nb;
            this.created = false(1,nb);
            hlp = reshape(repmat(1:nb,bcols,1),[],1);
            hlp2 = reshape(repmat(1:bcols,nb,1)',[],1);
            this.idx = [hlp(1:m) hlp2(1:m)];
            
            if nb == 1
                this.loadBlock(1);
            end
            % Matrix case: Assign value directly
            if matrixin
                this.subsasgn(struct('type',{'()'},'subs',{{':',':'}}),A);
            end
        end
        
        function B = spawnWithContent(this, A)
            % Creates a new FileMatrix containing the matrix A. The matrix is stored at the 
            % at the same location as the current matrix and the same block size is used.
            %
            % Parameters:
            % A: The new matrix. @type matrix<double>
            %
            % Return values:
            % B: A new data.FileMatrix instance containing B.
            B = data.FileMatrix(A,'Dir',fileparts(this.DataDirectory),'Blocksize',this.blocksize);
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
        
        function pos = getBlockPos(this, nr)
            % Returns the column indices of the block "nr" within the full matrix.
            pos = (nr-1)*this.bCols + 1 : min(nr*this.bCols,this.m);
        end
        
        function copy = copyWithNewBlockSize(this, block_size)
            copy = data.FileMatrix(this.n,this.m,...
                'Dir',fileparts(this.DataDirectory),'BlockSize',block_size);
            for k=1:this.nBlocks
                pos = this.getBlockPos(k);
                copy.subsasgn(struct('type',{'()'},'subs',{{':',pos}}),this.loadBlock(k));
            end
        end
        
        function relocate(this, new_root)
            % Relocates this FileMatrix instance to a new folder
            %
            % Note that the FileMatrix takes a directory 'inside' which it will create a new
            % directory using a certain hash.
            % Thus, upon relocation, the folder that 'contains' the created folder must be
            % passed.
            %
            % For example, assume we have DataDirectory=/some/path/matrix_246sg351dg. Then
            % relocate with '/some/new/root' leads to the new path
            % '/some/new/root/matrix_246sg351dg' at which the blocks of this file matrix are
            % assumed to reside.
            %
            % Parameters:
            % new_root: The root folder that contains a FileMatrix self-defined folder used to
            % store the matrix blocks.
            %
            % See also: data.FileData.relocate
            %
            % @new{0,7,dw,2013-05-28} Added this method.
            
            % Extract matrix_xxx folder name
            [~, mfolder] = fileparts(this.DataDirectory);
            % Concatenate with new path and call superclass relocate
            relocate@data.FileData(this, fullfile(new_root, mfolder));
        end
        
        function res = transposedTimes(this, B)
            % Implements the operation A'*B for this matrix A and another FileMatrix B
            if ~isa(B,'data.FileMatrix')
                error('Operation only applicable to other FileMatrix instances');
            end
            res = data.FileMatrix(this.m,B.m,'Dir', fileparts(B.DataDirectory),...
                    'BlockSize', B.blocksize);
            key = struct('type',{'()'},'subs',[]);
            for j = 1:B.nBlocks
                Bb = B.getBlock(j);
                for i=1:this.nBlocks
                    key.subs = {this.getBlockPos(i),B.getBlockPos(j)};
                    res.subsasgn(key,this.getBlock(i)'*Bb);
                end
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
            % Special case for one block
            if this.nBlocks == 1
                this.ensureCacheLoaded;
                value = this.cachedBlock.^expo;
                return;
            end
            value = data.FileMatrix(this.n,this.m,'Dir',fileparts(this.DataDirectory),...
                'BlockSize', this.blocksize);
            for k=1:this.nBlocks
                value.saveBlock(k,this.loadBlock(k).^expo);
            end
        end
        
        function n = numel(~)
            n = 1;
        end
        
        function [n, m] = size(this, dim)
            % Implementation of ABlockedData.size
            n = [this.n this.m];
            if nargin == 2
                if dim > 0 && dim < 3
                    n = n(dim);
                else
                    n = 0;
                end
            elseif nargout == 2
                n = this.n;
                m = this.m;
            end
        end
        
        function varargout = subsref(this, key)
            % Implements subscripted value retrieval.
            %
            % See also: subsref
            if strcmp(key(1).type,'()')
                if this.nBlocks == 1
                    this.ensureCacheLoaded;
                    [varargout{1:nargout}] = builtin('subsref', this.cachedBlock, key);
                    return;
                end
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
                    value = zeros(the_n,size(pos,1));
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
                if this.nBlocks == 1
                    this.ensureCacheLoaded;
                    this.cachedBlock = builtin('subsasgn', this.cachedBlock, key, value);
                    this.cacheDirty = true;
                    return;
                end
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
                this = builtin('subsasgn', this, key, value);
            end
        end
        
        function diff = minus(A, B)
            if all(size(A) == size(B))
                if isa(A,'data.FileMatrix')
                    if isa(B,'data.FileMatrix')
                        error('Not yet implemented.');
                        %                         diff = data.FileMatrix(A.n, A.m, A);
                        %                         key = struct('type',{'()'},'subs',{{':'}});
                        %                         for i=1:A.nBlocks
                        %                             pos = A.getBlockPos(i);
                        %                             key.subs{2} = pos;
                        %                             diff.subsasgn(key,A.loadBlock(i) - B(:,pos));
                        %                         end
                    else
                        % Special case for one block
                        if A.nBlocks == 1
                            A.ensureCacheLoaded;
                            diff = A.cachedBlock - B;
                            return;
                        end
                        %| @todo need FileMatrix as result type here
                        diff = zeros(A.n, A.m);
                        for i=1:A.nBlocks
                            pos = A.getBlockPos(i);
                            diff(:,pos) = A.loadBlock(i) - B(:,pos);
                        end
                    end
                else
                    diff = minus(B, A);
                end
            else
                error('Invalid matrix dimensions.');
            end
        end
        
        function a = mldivide(L, R)
            if isa(L,'data.FileMatrix')
                if L.nBlocks == 1
                    L.ensureCacheLoaded;
                    a = L.cachedBlock\R;
                else
                    error('Not yet implemented');
                end
            else
                if R.nBlocks == 1
                    R.ensureCacheLoaded;
                    a = L\R.cachedBlock;
                else
                    error('Not yet implemented');
                end
            end
        end
        
        function AB = mtimes(A, B)
            % Override of ABlockedData.mtimes
            if isa(A,'data.FileMatrix')
                % FileMatrix * FileMatrix case
                if isa(B,'data.FileMatrix')
                    if A.nBlocks == 1 && B.nBlocks == 1
                        A.ensureCacheLoaded;
                        B.ensureCacheLoaded;
                        AB = data.FileMatrix(A.cachedBlock * B.cachedBlock);
                        return;
                    else
                        % Set blocksize such that the result has the same number of columns per
                        % block as B has
                        AB = data.FileMatrix(A.n,B.m,'Dir', fileparts(B.DataDirectory),...
                        'BlockSize', A.blockSizeOf(A.n,B.bCols));
                        key = struct('type',{'()'},'subs',[]);
                        for i=1:B.nBlocks
                            Bb = B.getBlock(i);
                            key.subs = {':',B.getBlockPos(i)};
                            for k = 1:A.nBlocks
                                Ab = A.getBlock(k);
                                sum = AB.subsref(key) + Ab*Bb(A.getBlockPos(k),:);
                                AB.subsasgn(key,sum);
                            end
                        end
                    end
                else
                    % Special case for one block
                    if A.nBlocks == 1
                        A.ensureCacheLoaded;
                        AB = A.cachedBlock * B;
                        return;
                    end
                    if isscalar(B)  % Matrix*scalar case
                        AB = mtimes(B,A);
                    else
                        if A.m ~= size(B,1)
                            error('Matrix dimensions must agree.');
                        end
                        % FileMatrix * vec case
                        if size(B,2) == 1
                            AB = zeros(A.n,1);
                            for bidx = 1:A.nBlocks
                                % Only multiply if nonzero
                                if A.created(bidx)
                                    ABlock = A.loadBlock(bidx);
                                    AB = AB + ABlock*B(A.getBlockPos(bidx),:);
                                end
                            end
                            % FileMatrix * matrix case
                        else
                            AB = data.FileMatrix(A.n, size(B,2), 'Dir', fileparts(A.DataDirectory),...
                                'BlockSize', A.blocksize);
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
                end
                % [vec|mat] * FileMatrix case
            elseif isa(B,'data.FileMatrix')
                % Special case for one block
                if B.nBlocks == 1
                    B.ensureCacheLoaded;
                    AB = A*B.cachedBlock;
                    return;
                end
                if isscalar(A)  % scalar*Matrix
                    AB = data.FileMatrix(B.n, B.m, 'Dir', fileparts(B.DataDirectory),...
                        'BlockSize', B.blocksize);
                    key = struct('type',{'()'},'subs',{{':'}});
                    for bidx = 1:B.nBlocks
                        % Only multiply if nonzero
                        if B.created(bidx)
                            key.subs{2} = B.getBlockPos(bidx);
                            AB.subsasgn(key,A*B.loadBlock(bidx));
                        end
                    end
                else
                    AB = data.FileMatrix(size(A,1), B.m, 'Dir', fileparts(B.DataDirectory),...
                        'BlockSize', B.blocksize);
                    key = struct('type',{'()'},'subs',{{':'}});
                    % If AB has more blocks than B, the resulting matrix is larger than B.
                    % Thus, we need A*b to be as most as big as blocksize, which is why slices of B
                    % in the size of ABs blocks are taken.
                    if AB.nBlocks > B.nBlocks
                        for i=1:AB.nBlocks
                            key.subs{2} = AB.getBlockPos(i);
                            b = B.subsref(key);
                            AB.subsasgn(key,A*b);
                        end
                    else
                        for i=1:B.nBlocks
                            if B.created(i)
                                key.subs{2} = B.getBlockPos(i);
                                AB.subsasgn(key,A*B.loadBlock(i));
                            end
                        end
                        % Return a matrix if the A argument was a transposed vector
                        if size(A,1) == 1
                            AB = AB.toMemoryMatrix;
                        end
                    end
                end
            end
        end
        
        function AB = times(A, B)
            % Override of ABlockedData.mtimes
            if isa(A,'data.FileMatrix')
                % FileMatrix * FileMatrix case
                if isa(B,'data.FileMatrix')
                    error('Not yet implemented.');
                else
                    % Special case for one block
                    if A.nBlocks == 1
                        A.ensureCacheLoaded;
                        AB = A.cachedBlock .* B;
                        return;
                    elseif isscalar(B)
                        % Call matrix multiplication with scalar value
                        AB = mtimes(B,A);
                    else
                        if  A.n ~= size(B,1) || A.m ~= size(B,2)
                            error('Matrix dimensions must agree.');
                        end
                        AB = data.FileMatrix(A.n,A.m,'Dir', fileparts(A.DataDirectory),...
                            'BlockSize', A.blocksize);

                        for k=1:A.nBlocks
                            AB.saveBlock(k,A.loadBlock(k).*B(:,A.getBlockPos(k)));
                        end
                    end
                end
            else
                % Pointwise operation does not depend on which side it is done from.
                AB = times(B,A);
            end
        end
        
        function trans = ctranspose(this)
            trans = data.FileMatrix(this.m, this.n, 'BlockSize', this.blocksize, ...
                'Dir', fileparts(this.DataDirectory));
            % Special case for only one block
            if this.nBlocks == 1
                this.ensureCacheLoaded;
                trans = this.cachedBlock';
                return;
            end
            % Else: multiple blocks
            if this.InPlaceTranspose
                if this.nBlocks > 1
                    pi = ProcessIndicator('Creating in-place transposed of %d-block matrix in %s' ,this.nBlocks,...
                        false,this.nBlocks,trans.DataDirectory);
                end
                key = struct('type',{'()'},'subs',{{[],':'}});
                for j=1:this.nBlocks
                    B = this.loadBlock(j);
                    key.subs{1} = this.getBlockPos(j);
                    trans.subsasgn(key,B');
                    if this.nBlocks > 1
                        pi.step;
                    end
                end
            else
                if this.nBlocks > 1
                    pi = ProcessIndicator('Creating transposed of %d-block matrix in %s',2*this.nBlocks,...
                        false,this.nBlocks,trans.DataDirectory);
                end
                % Write out blocks in small chunks
                for j=1:this.nBlocks
                    B = this.loadBlock(j);
                    for k=1:trans.nBlocks
                        pos = trans.getBlockPos(k);
                        chunk = B(pos,:); %#ok
                        save(fullfile(trans.DataDirectory,sprintf('tmp_%d_%d.mat',j,k)), 'chunk');
                    end
                    if this.nBlocks > 1
                        pi.step;
                    end
                end
                % Read in blocks from respective chunks
                for k=1:trans.nBlocks
                    B = trans.loadBlock(k);
                    for j=1:this.nBlocks
                        file = fullfile(trans.DataDirectory,sprintf('tmp_%d_%d.mat',j,k));
                        s = load(file);
                        pos = this.getBlockPos(j);
                        B(pos,:) = s.chunk';
                        % remove temp file
                        delete(file);
                    end
                    trans.saveBlock(k,B);
                    if this.nBlocks > 1
                        pi.step;
                    end
                end
            end
            if this.nBlocks > 1
                pi.stop;
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
                        if A.created(i) || A.cachedNr == i
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
            if ~isempty(this.created) && ~this.isSaved
                % Remove exactly the files for the matrix
                for k=1:this.nBlocks
                    if this.created(k)
                        f = fullfile(this.DataDirectory, sprintf('block_%d.mat',k));
                        if exist(f,'file')
                            delete(f);
                        end
                    end
                end
            end
            % Superclass delete removes the folder if empty.
            delete@data.FileData(this);
        end
        
        %% data.ABlockedData implementations
        function n = getNumBlocks(this)
            % Implementation of ABlockedData.getNumBlocks
            n = this.nBlocks;
        end
        
        function B = getBlock(this, nr)
            % Implementation of ABlockedData.getBlock
            B = this.loadBlock(nr);
        end
    end
    
    methods(Access=protected)
        
        function this = saveobj(this)
            this = saveobj@data.FileData(this);
            if this.cacheDirty
                A = this.cachedBlock;%#ok
                save([this.DataDirectory filesep sprintf('block_%d.mat',this.cachedNr)], 'A');
                this.created(this.cachedNr) = true;
            end
            this.cachedBlock = [];
            this.cachedNr = [];
            this.cacheDirty = false;
        end
    end
    
    methods(Access=private)
        
        function ensureCacheLoaded(this)
            if this.nBlocks == 1 && isempty(this.cachedBlock)
                this.loadBlock(1);
            end
        end
        
        function saveBlock(this, nr, A)
            % Also save changes to cachedBlock if number matches
            if nr == this.cachedNr
                this.cachedBlock = A;
                this.cacheDirty = true;
            else
                % Save actual block that should be saved
                save([this.DataDirectory filesep sprintf('block_%d.mat',nr)], 'A');
            end
            this.updateMinMax(A);
            this.created(nr) = true;
        end
        
        function A = loadBlock(this, nr)
            if nr == this.cachedNr
                A = this.cachedBlock;
            else
                % Before loading a new block, save the old one if the cache is dirty!
                if this.cacheDirty
                    A = this.cachedBlock;
                    save([this.DataDirectory filesep sprintf('block_%d.mat',this.cachedNr)], 'A');
                    this.updateMinMax(A);
                    this.created(this.cachedNr) = true;
                    % cacheDirty is set "false" at end!
                end
                if this.created(nr)
                    s = load([this.DataDirectory filesep sprintf('block_%d.mat',nr)]);
                    A = s.A;
                else
                    % Ensures correct size for block (even if only one with less than bCols
                    % columns)
                    A = zeros(this.n, length(this.getBlockPos(nr)));
                end
                this.cachedNr = nr;
                this.cachedBlock = A;
                this.cacheDirty = false;
            end
        end
        
        function updateMinMax(this,A)
            this.MinValue = min(this.MinValue,min(min(A)));
            this.MaxValue = max(this.MaxValue,max(max(A)));
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this, initfrom)
            % Loads a FileMatrix instance.
            
            created = false;
            if ~isa(this, 'data.FileMatrix')
                initfrom = this;
                this = data.FileMatrix(initfrom.n,initfrom.m,'Dir',...
                    initfrom.DataDirectory,'BlockSize',initfrom.blocksize);
                created = true;
            end
            if nargin == 2 || created
                this.InPlaceTranspose = initfrom.InPlaceTranspose;
                this.bCols = initfrom.bCols;
                this.nBlocks = initfrom.nBlocks;
                this.n = initfrom.n;
                this.m = initfrom.m;
                this.MinValue = initfrom.MinValue;
                this.MaxValue = initfrom.MaxValue;
                this.idx = initfrom.idx;
                this.created = initfrom.created;
                this.blocksize = initfrom.blocksize;
                if this.nBlocks == 1
                    field = 'cachedBlock';
                    if isfield(initfrom,'block1')
                        field = 'block1';
                    end
                    this.cachedBlock = initfrom.(field);
                end
                this = loadobj@data.FileData(this, initfrom);
            else
                this = loadobj@data.FileData(this);
            end
        end
    end
    
    methods(Static)
        function fm = recoverFrom(directory)
            % Tries to recover a FileMatrix from a given directory, containing the old block
            % data files.
            %
            % Parameters:
            % directory: The directory to recover from. If not given, KerMor.getDir is used.
            % @type char @default KerMor.getDir
            %
            % Return values:
            % fm: The recovered FileMatrix instance. [] if the KerMor.getDir dialog has been
            % aborted. @type data.FileMatrix
            if nargin == 0
                directory = KerMor.getDir;
                if ~directory
                    fm = [];
                    return;
                end
            end
            d = dir(directory);
            nf = length(d)-2;
            if nf > 0
                nb = 0;
                for k=1:nf
                    if ~exist(fullfile(directory,sprintf('block_%d.mat',k)),'file')
                        break;
                    end
                    nb = nb + 1;
                end
                if nb > 0
                    s = load(fullfile(directory,'block_1.mat'));
                    A1 = s.A;
                    [N, nCols] = size(A1);
                    if nb > 1
                        s = load(fullfile(directory,sprintf('block_%d.mat',nb)));
                    end
                    M = (nb-1)*nCols + size(s.A,2);
                    fm = data.FileMatrix(N,M,'Dir',fileparts(directory),'BlockSize',8*nCols*N);
                    fm.updateMinMax(A1); % first block
                    if nb > 1
                        fm.updateMinMax(s.A);
                        for k=2:nb-1
                            s2 = load(fullfile(directory,sprintf('block_%d.mat',k)));
                            fm.updateMinMax(s2.A);
                        end
                    end
                    % Delete & reassign correct directory
                    rmdir(fm.DataDirectory);
                    fm.DataDirectory = directory;
                    fm.created(:) = true;
                    fm.isSaved = true;
                else
                    error('No block_%%.mat files of a FileMatrix found.');
                end
            else
                error('No files found in %s',directory);
            end
        end
        
        function bs = blockSizeOf(arg1, arg2)
            % Computes the block size in Megabytes for a given matrix of matrix
            % dimensions.
            %
            % The matrices are assumed to contain double values.
            %
            % Parameters:
            % arg1: Either a double matrix or the row number
            % arg2: The column number @type integer @default []
            if nargin < 2
                n = numel(arg1);
            else
                n = arg1*arg2;
            end
            bs = n*8/1024^2;
        end
        
    end
end