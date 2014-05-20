classdef Norm
% Norm: Static class for commonly used norms on sets of vectors
%
% All norm functions perform their respective computations on each column
% of a given matrix!
%
% @author Daniel Wirtz @date 2012-03-29
%
% @new{0,6,dw,2012-03-29} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods(Static)
        function n = L2(x)
            % Returns the discrete `L^2` norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix `\vX` containing column vectors `\vx_i` @type matrix<double>
            %
            % Return values:
            % n: The `L^2` norm `\norm{\vx}{2} = \sqrt{\sum\limits_{i=1}^dx_i^2}` for each column
            % vector `\vx\in\R^d` in `\vX` @type rowvec<double>
            n = sqrt(sum(x.^2,1));
        end
        
        function vecs = normalizeL2(vecs)
            no = Norm.L2(vecs);
            vecs = vecs ./ (ones(1,size(vecs,1))' * no);
        end
        
        function n = L1(x)
            % Returns the discrete `L^1` norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix `\vX` containing column vectors `\vx_i` @type matrix<double>
            %
            % Return values:
            % n: The `L^1` norm `\norm{\vx}{1} = \sum\limits_{i=1}^d|x_i|` for each column vector
            % `\vx\in\R^d` in `\vX` @type rowvec<double>
            n = sum(abs(x),1);
        end
        
        function n = Linf(x)
            % Returns the discrete `L^\infty` norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix `\vX` containing column vectors `\vx_i` @type matrix<double>
            %
            % Return values:
            % n: The `L^\infty` norm `\norm{\vx}{\infty} = \max\limits_{i=1}^d|x_i|` for each
            % column vector `\vx\in\R^d` in `\vX` @type rowvec<double>
            n = max(abs(x),[],1);
        end
        
        function n = LG(x, G)
            % Returns the `\vG`-induced norm for each column vector in `\vX`.
            %
            % Parameters:
            % x: A matrix `\vX` containing column vectors `\vx_i` @type matrix<double>
            % G: A positive definite matrix `\vG`
            %
            % Return values:
            % n: The G-norm `\noG{\vx} = \sqrt{\vx^T\vG\vx}` for each column vector `\vx` in
            % `\vX` @type rowvec<double>
            n = sqrt(sum(x.*(G*x),1));
        end
    end
end