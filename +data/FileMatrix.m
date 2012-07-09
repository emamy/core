classdef FileMatrix < data.FileData
% FileMatrix: 
%
%
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
        
        % The machine format (little/big endian, 32/64bit) to use for binary files.
        %
        % @default ieee-be (Big Endian)
        %
        % See 'doc fopen' for more details.
%         MachineFormat = 'ieee-be';
    end

    properties%(Access=private)
        n,m;
        
        bRows;
        
        nBlocks;
        
        idx;
        
        created;
    end
    
    methods
        function this = FileMatrix(n, m, storage_root, block_size)
            % Parameters:
            % block_size: The maximum block size in Bytes. @default 256MB
            if nargin < 4
                block_size = min(data.FileMatrix.BLOCK_SIZE,n*m*8);
                if nargin < 3
                   storage_root = KerMor.App.DataStoreDirectory;
                end
            end
            ID = general.IDGenerator.generateID;
            this = this@data.FileData(fullfile(storage_root,['matrix_' num2str(ID)]));
            this.n = n; 
            this.m = m;
            this.bRows = floor(block_size/(8*m));
            this.nBlocks = ceil(n / this.bRows);
            this.created = false(1,this.nBlocks);
            hlp = reshape(repmat(1:this.nBlocks,this.bRows,1),[],1);
            hlp2 = reshape(repmat(1:this.bRows,this.nBlocks,1)',[],1);
            this.idx = [hlp(1:n) hlp2(1:n)];
        end
        
        function [V,S,W] = getSVD(this, k)
            psize = min(this.n,this.m);
            if nargin < 2
                k = psize;
            end
            
            opts.issym = 1;
            colmode = this.n < this.m;
            fprintf('FileMatrix: Computing truncated %d-SVD on %dx%d matrix (%d blocks)...\n',...
                    k,this.n,this.m,this.nBlocks);
            if colmode
                [V,S] = eigs(@colmode_mult,psize,k,'la',opts);
                if nargout > 2
                    W = (this')*(V/S);
                end
            else
                [W,S] = eigs(@rowmode_mult,psize,k,'la',opts);
                V = this*(W/S);
            end
            
            function w = rowmode_mult(v)
                w = 0;
                for j = 1:this.nBlocks
                    B = this.loadBlock(j);
                    w = w + B'*(B*v);
                end
            end
            
            function w = colmode_mult(v)
                w = zeros(size(v));
                for j = 1:this.nBlocks
                    B = this.loadBlock(j);
                    pos = (j-1)*this.bRows + 1 : min(j*this.bRows,size(v,1));
                    w(pos,:) = B*(B'*v(pos,:));
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
                    if strcmp(s{1},':')
                        s{1} = 1:this.n;
                    end
                    pos = this.idx(s{1},:);
                    blocks = unique(pos(:,1));
                    m = this.m;
                    if ~strcmp(s{2},':')
                        m = length(s{2});
                    end
                    value = zeros(length(s{1}),m);
                    for bidx = 1:length(blocks)
                        b = blocks(bidx);
                        B = this.loadBlock(b);
                        value(pos(:,1) == b,:) = B(pos(pos(:,1)==b,2),s{2});
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
                    if strcmp(s{1},':')
                        s{1} = 1:this.n;
                    end
                    pos = this.idx(s{1},:);
                    blocks = unique(pos(:,1));
                    for bidx = 1:length(blocks)
                        b = blocks(bidx);
                        B = this.loadBlock(b);
                        B(pos(pos(:,1)==b,2),s{2}) = value(pos(:,1) == b,:);
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
                        pos = (bidx-1)*this.bRows + 1 : min(bidx*this.bRows,this.n);
                        Ax(pos,:) = B*x;
                    end
                end
            end
        end
        
        function trans = ctranspose(this)
            root = fileparts(this.DataDirectory);
            bs = 8*this.bRows*this.m;
            trans = data.FileMatrix(this.m,this.n,root,bs);
            for j=1:this.nBlocks
                B = this.loadBlock(j);
                pos = (j-1)*this.bRows + 1 : min(j*this.bRows,this.n);
                % Why EVER this doesnt work.. GRRR
                %trans(:,pos) = B';
                trans.subsasgn(struct('type',{'()'},'subs',{{':', pos}}),B');
            end
        end
    end
    
    methods(Access=private)
        
        function saveBlock(this, nr, A)%#ok
            file = fullfile(this.DataDirectory, ['block_' num2str(nr) '.mat']);
            save(file, 'A');
            this.created(nr) = true;
        end
        
        function A = loadBlock(this, nr)
            if this.created(nr)
                file = this.getfile(['block_' num2str(nr) '.mat']);
                s = load(file);
                A = s.A;
            else
                A = zeros(this.bRows,this.m); 
            end
        end
    end
end