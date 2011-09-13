classdef MemoryKernelMatrix < data.AKernelMatrix
% MemoryKernelMatrix: A simple implementation of a AKernelMatrix wrapping
% around a memory-based kernel matrix and forwarding calls.
%
% @todo Extra layer seems expensive in terms of computational cost,
% consider removing it.. but it's only offline phase.
%
% @author Daniel Wirtz @date 2011-09-09
%
% @new{0,5,dw,2011-09-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetObservable, Dependent)
        % Flag whether to build the inverse matrix upon augmentation of the
        % kernel matrix.
        %
        % @propclass{optional} Optional property that can be used in
        % algorithms.
        %
        % @default false
        BuildInverse;
       
        % Flag that indicates whether to keep an LU decomposition of the kernel
        % matrix or not.
        %
        % @propclass{optional} Speedup maybe gained if subsequent calls to interpolation using the
        % same kernel matrix are made.
        %
        % @default false
        UseLU;
    end
    
    properties(SetAccess=private)
        % The inner kernel matrix
        K;
        
        % The inverse matrix
        Kinv = [];
        
        % Private variable to store the lower left part of the optional LU
        % decomp.
        %
        % See also: K UseLU
        L;
        
        % Private variable to store the lower left part of the optional LU
        % decomp.
        %
        % See also: K UseLU
        U;
    end
    
    properties(Access=private)
        fBInv = false;
        fUseLU = false;
    end
    
    methods
        function this = MemoryKernelMatrix(K)
            if nargin == 1
                this.K = K;
            end
        end
        
        function this = augment(this, vec)
            % @todo performance check - see if preallocating and direct
            % assignment is faster than two matrix compositions
            %
            % @todo check if there is a way to update the LU instead of
            % recomputation?
            vec = reshape(vec,[],1);
            this.K = [[this.K; vec(1:end-1)'] vec];
            if this.fUseLU
                [this.L, this.U] = lu(this.K);
            end
            if this.fBInv
                if numel(this.K) == 1
                    this.Kinv = 1/this.K;
                else
                    this.Kinv = general.MatUtils.getExtendedInverse(this.Kinv, vec);
                end
            end
        end
        
        function value = subsref(this, key)
            % Implements subscripted value retrieval.
            %
            % See also: subsref
            if strcmp(key(1).type,'.')
                value = builtin('subsref', this, key);
            else
                value = builtin('subsref', this.K, key);
            end
            
        end
        
        function this = subsasgn(this, key, value)
            % Implements subscripted assignment.
            %
            % See also: subsasgn
            if strcmp(key(1).type,'.')
                builtin('subsasgn', this, key, value);
            else
                builtin('subsasgn', this.K, key, value);
            end
        end
        
        function prod = mtimes(A, B)
            % Implements the default multiplication method.
            if isa(A,'data.MemoryKernelMatrix') && isa(B,'data.MemoryKernelMatrix')
                prod = data.MemoryKernelMatrix(A.K * B.K);
            elseif isa(A,'data.MemoryKernelMatrix')
                if isscalar(B)
                    prod = data.MemoryKernelMatrix(A.K * B);
                else
                    prod = A.K * B;
                end
            elseif isa(B,'data.MemoryKernelMatrix')
                if isscalar(A)
                    prod = data.MemoryKernelMatrix(A * B.K);
                else
                    prod = A * B.K;
                end
            end
        end
        
        function a = mldivide(this, vec)
            a = this.K\vec;
        end
        
        function value = size(this)
            value = size(this.K);
        end
        
        function set.BuildInverse(this, value)
            this.fBInv = value;
            if isempty(this.Kinv) && ~isempty(this.K)
                this.Kinv = inv(this.K);
            end
        end
        
        function value = get.BuildInverse(this)
            value = this.fBInv;
        end
        
        function set.UseLU(this, value)
            % Sets the UseLU property and ensures that if it was set after
            % the kernel matrix was assigned the LU decomposition is
            % computed in each case.
            this.fUseLU = value;
            if value && ~isempty(this.K)
                [this.L, this.U] = lu(this.K);
            end
        end
        
        function flag = get.UseLU(this)
            flag = this.fUseLU;
        end
    end
    
end