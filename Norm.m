classdef Norm
% Norm: Static class for commonly used norms
%
% @author Daniel Wirtz @date 2012-03-29
%
% @new{0,6,dw,2012-03-29} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        function n = L2(x)
            % Returns the discrete `L^2` norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix containing column vectors @type matrix<double>
            %
            % Return values:
            % n: The `L^2` norm `||v||_2 = \sqrt{\sum\limits_{i=1}^d}v_i^2}` for each column
            % vector `v\in\R^d` @type rowvec<double>
            n = sqrt(sum(x.^2,1));
        end
        
        function n = L1(x)
            % Returns the discrete `L^1` norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix containing column vectors @type matrix<double>
            %
            % Return values:
            % n: The `L^1` norm `||v||_1 = \sum\limits_{i=1}^d|v_i|` for each column vector
            % `v\in\R^d` @type rowvec<double>
            n = sum(abs(x),1);
        end
        
        function n = Linf(x)
            % Returns the discrete `L^\infty` norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix containing column vectors @type matrix<double>
            %
            % Return values:
            % n: The `L^\infty` norm `||v||_\infty = \max\limits_{i=1}^d|v_i|` for each column
            % vector `v\in\R^d` @type rowvec<double>
            n = max(abs(x),[],1);
        end
        
        function n = LG(x, G)
            % Returns the G-induced norm for each column vector in x.
            %
            % Parameters:
            % x: A matrix containing column vectors @type matrix<double>
            % G: A positive definite matrix `G`
            %
            % Return values:
            % n: The G-norm `||v||_G = \sqrt{x^tGx}` for each column vector `v` @type rowvec<double>
            n = sqrt(sum(x.*(G*x),1));
        end
    end
    
end