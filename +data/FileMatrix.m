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
    end
    
    properties
        % The minimum relative value of singular values that triggers selection of the 
        % compared to the largest one.
        MinRelSingularValueSize = 1e-12;
    end

    properties(SetAccess=private)
        n,m;
        
        bRows;
        
        nBlocks;
    end
    
    properties(Access=private)
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
            this.bRows = max(floor(block_size/(8*m)),1);
            this.nBlocks = ceil(n / this.bRows);
            this.created = false(1,this.nBlocks);
            hlp = reshape(repmat(1:this.nBlocks,this.bRows,1),[],1);
            hlp2 = reshape(repmat(1:this.bRows,this.nBlocks,1)',[],1);
            this.idx = [hlp(1:n) hlp2(1:n)];
        end
        
        function [Vs,S,Vl] = getSVD(this, k)
            psize = min(this.n,this.m);
            if nargin < 2
                k = psize;
            end
            
            opts.issym = 1;
            colmode = this.n < this.m;
            fprintf('FileMatrix: Computing truncated %d-SVD on %dx%d matrix (%d blocks)...\n',...
                    k,this.n,this.m,this.nBlocks);
            
            fun = @rowmode_mult;
            if colmode
                warning('KerMor:FileMatrix',['Computing SVD on matrix with n < m is very inefficient. '...
                    'Consider arranging the matrix transposed.']);
                fun = @colmode_mult;    
            end
            [Vs,S] = eigs(fun,psize,k,'la',opts);
            sel = sqrt(diag(S)/S(1)) >= this.MinRelSingularValueSize;
            Vs = Vs(:,sel);
            if size(Vs,2) < k
                warning('KerMor:FileMatrix','Have only %d nonzero singular values instead of %d desired ones.',...
                    size(Vs,2),k);
            end
            S = sqrt(S(sel,sel));
            if nargout > 2
                if colmode
                    Vl = this'*(Vs/S);
                else
                    Vl = this*(Vs/S);
                end
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
                    Bj = this.loadBlock(j);
                    for i = 1:this.nBlocks
                        Bi = this.loadBlock(i);
                        posi = (i-1)*this.bRows + 1 : min(i*this.bRows,size(v,1));
                        posj = (j-1)*this.bRows + 1 : min(j*this.bRows,size(v,1));
                        w(posj,:) = w(posj,:) + Bj*(Bi'*v(posi,:));
                    end
                end
            end
        end
        
        function n = numel(~)
            n = 1;
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
                trans.subsasgn(struct('type',{'()'},'subs',{{':', pos}}),B');
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
    end
    
    methods(Access=private)
        function saveBlock(this, nr, A)%#ok
            save([this.DataDirectory filesep sprintf('block_%d.mat',nr)], 'A');
            this.created(nr) = true;
        end
        
        function A = loadBlock(this, nr)
            if this.created(nr)
                s = load([this.DataDirectory filesep sprintf('block_%d.mat',nr)]);
                A = s.A;
            else
                A = zeros(this.bRows,this.m);
                % Shorten the last block to effectively used size
                if nr == this.nBlocks
                    A = A(1:(this.n-(this.nBlocks-1)*this.bRows),:);
                end
            end
        end
    end
    
    methods(Static)
        function res = test_FileMatrix
            res = true;
            B = rand(10,100);
            A = data.FileMatrix(10,100,KerMor.App.DataStoreDirectory,1600);
            key = struct('type',{'()'},'subs',[]);
            for k=1:10
                key.subs = {k, ':'};
                A.subsasgn(key,B(k,:));
            end
            key.subs = {':', ':'};
            res = res && all(all(A.subsref(key) == B));
            
            % Transpose test
            At = A';
            res = res && all(all(At.subsref(key) == B'));
            
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
            [Vt,St,Ut] = At.getSVD;
            [Vt5,St5,Ut5] = At.getSVD(5);
            res = res && norm(abs(Vt)-abs(vt),'fro') < p && norm(abs(Ut)-abs(ut),'fro') < p &&...
                norm(diag(St)-diag(st)) < p;
            res = res && norm(abs(Vt5)-abs(vt(:,1:5)),'fro') < p ...
                    && norm(abs(Ut5)-abs(ut(:,1:5)),'fro') < p ...
                    && norm(diag(St5)-diag(st(1:5,1:5))) < p;
        end
    end
end