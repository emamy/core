classdef AffParamMatrix < ICloneable
    % General time/parameter-affine matrix
    %
    % Represents a linear combination of matrices `A_i` with scalar parameter and/or time dependent
    % coefficient functions `\theta_i` of the form `\sum\limits_{i=0}^Q \theta_i(t,\mu)A_i`.
    %
    % Some default operators (minus, mtimes, ctranspose) have been implemented to work on instances
    % of this type. This enables very easy notation and increases readability when used.
    %
    % @ingroup general
    %
    % @change{0,5,dw,2011-07-05} Renamed this class to AffParamMatrix and removed the IProjectable
    % interfaces. Renamed the 'evaluate' function to 'compose'.
    %
    % @change{0,4,sa,2011-05-06} Implemented Setters for the class
    % properties
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing    
    %
    % @todo implement/override other arithmetic operations
        
    properties(SetAccess=private)
        % The number of affine matrices / size of the linear combination
        N = 0;
        
        % The matrices for the affine function
        % @propclass{critical}
        Matrices = [];
    end
    
    properties(Access=private)
        funStr = {};
        cfun = [];
        dims;
    end
    
    methods 
        
        function M = compose(this, t, mu)
            % Composes the affine-linear combination of matrices for given time `t` and parameter
            % `\mu`.
            %
            % Parameters:
            % t: The current time `t`
            % mu: The current parameter vector `\mu`
            %
            % Return values:
            % M: The sum `M(t,\mu) = \sum\limits_{i=0}^Q\theta_i(t,\mu)A_i`
            if this.N == 0
                M = [];
                return;
            end
            M = reshape(sum(this.Matrices * this.cfun(t,mu),2),this.dims);
        end
        
%         function varargout = subsref(this, S)
%             if S(1).type(1) == '('
%                 varargout{1} = this.compose(S.subs{1},S.subs{2});
%             else
%                 [varargout{1:nargout}] = builtin('subsref',this,S);
%             end
%         end
                
        function copy = clone(this, copy)
            % Creates a copy of this affine parametric matrix.
            if nargin == 1
                copy = general.AffParamMatrix;
            end
            copy.Matrices = this.Matrices;
            copy.N = this.N;
            copy.funStr = this.funStr;
            copy.dims = this.dims;
            copy.buildCoeffFun;
        end
        
        function addMatrix(this, coeff_fun, mat)
            % Adds a matrix with corresponding coefficient function to the affine parametric matrix.
            %
            % Parameters:
            % coeff_fun: A string describing the `i`-th coefficient function's evaluation as if
            % entered into a function handle interior or a function. By convention, the passed arguments are named
            % 't' and 'mu'.
            % For example, one would pass '(t+1)/mu(1))' as 'coeff_fun' string if the function was to be
            % '@(t,mu)(t+1)/mu(1))'
            % mat: The corresponding matrix `A_i`
          
            if ~isa(coeff_fun,'char')
                error('Coeff_fun must be a string.');
            elseif ~isreal(mat)
                error('Parameter "mat" must be a real matrix.');
            end
            if isempty(this.dims)
                this.dims = size(mat);
            end
            this.N = this.N + 1;
            
            this.funStr{end+1} = coeff_fun;
            %this.Coefficients{end+1} = coeff_fun;
            
            this.Matrices(:,end+1) = mat(:);
            
            % Create the coefficient function
            this.buildCoeffFun;
        end
        
        function pr = mtimes(A, B)
            % Implements the default multiplication method.
            if isa(A,'general.AffParamMatrix') && isa(B,'general.AffParamMatrix')
                pr = general.AffParamMatrix;
                pr.N = A.N * B.N;
                pr.dims = [A.dims(1) B.dims(2)];
                pr.Matrices = zeros(prod(pr.dims),pr.N);
                for i = 1:A.N
                    M = reshape(A.Matrices(:,i),A.dims);
                    for j = 1:B.N
                        pr.Matrices(:,(i-1)*B.N + j) = ...
                            reshape(M * reshape(B.Matrices(:,j),B.dims),1,[]);
                        pr.funStr{end+1} = ['(' A.funStr{i} ')*(' B.funStr{j} ')'];
                    end
                end
                pr.buildCoeffFun;
            elseif isa(A,'general.AffParamMatrix')
                pr = A.clone;
                for i=1:A.N
                    pr.Matrices(:,i) = reshape(...
                        reshape(A.Matrices(:,i),A.dims) * B,1,[]);
                end 
            elseif isa(B,'general.AffParamMatrix')
                pr = B.clone;
                for i=1:B.N
                    pr.Matrices(:,i) = reshape(...
                        A * reshape(B.Matrices(:,i),B.dims),1,[]);
                end 
            end
        end
        
        function diff = minus(A, B)
            % Implements the default substraction method.
            %
            % So far only works for the case of identical coefficient functions, as more
            % sophisticated differences for different coefficient functions and combination sizes
            % are not used yet.
            if isa(A,'general.AffParamMatrix') && isa(B,'general.AffParamMatrix')
                if A.N == B.N && all(strcmp(A.funStr,B.funStr))
                    diff = A.clone;
                    diff.Matrices = A.Matrices - B.Matrices;
                else
                    error('If two AffParamMatrices are subtracted, the number of elements must be the same and the coefficient functions must match.');
                end
            elseif isa(A,'general.AffParamMatrix')
                diff = A.clone;
                diff.addMatrix('-1',B)
            elseif isa(B,'general.AffParamMatrix')
                diff = B.clone;
                diff.Matrices = -diff.Matrices;
                diff.addMatrix('1',A);
            end
        end
        
        function transp = ctranspose(this)
            % Implements the transposition for affine parametric matrices.
            transp = this.clone;
            for i=1:this.N
                transp.Matrices(:,i) = reshape(reshape(this.Matrices(i),this.dims)',1,[]);
            end
        end
        
    end
    
    methods(Access=private)
        function buildCoeffFun(this)
            % Creates the coefficient function handle from the recorded funStr values. This function
            % was exported here in order to create a new function handle upon cloning, which ensures
            % a clean workspace for the new handle at cloning.
            this.cfun = eval(['@(t,mu)[' general.Utils.implode(this.funStr,'; ') ']']);
        end
    end
    
end

