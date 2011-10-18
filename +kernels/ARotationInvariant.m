classdef ARotationInvariant < KerMorObject
    % Abstract class for rotation invariant kernels
    %
    % All rotation invariant kernels have the form
    % `` \Phi(x,y) := \phi(||x-y||), x\in\mathbb{R}^d``
    % for some real-valued scalar function `\phi`.
    %
    % When combinations of Kernels are used, this interface will have to be changed to a property.
    % Up to now, the class CombinationKernel cannot dynamically adopt to the interface for the case
    % that all contained kernels implement this interface.
    %
    % @attention The evaluation of any rotation invariant kernel bases on
    % fast computation of `||x-y||^2`, albeit only `||x-y||` is passed to
    % the ARotationInvariant.evaluateScalar function (as it is expected). However, should the
    % ARotationInvariant.evaluateScalar function square the argument again, it should be
    % considered overriding the ARotationInvariant.evaluate function in this class to improve
    % performance. An example for this is the GaussKernel.
    %
    % @author Daniel Wirtz @date 2011-08-09
    %
    % @new{0,5,dw,2011-10-17} 
    % - Added this class.
    % - Implemented the general evaluate function for rotation invariant
    % kernels.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    %
    % @todo: change property to check for interface implementation,
    % implement in other suitable kernels
    
    properties(SetObservable)
        % The norm to use for the evaluation of `||x-y||_G`.
        %
        % Must be a positive definite, symmetric matrix `G`.
        %
        % @propclass{critical} If a custom norm is used (i.e. after
        % subspace projection) tihs must be set in order to obtain correct
        % evaluations.
        %
        % @type matrix @default 1
        G = 1;
    end
    
    methods
        function this = ARotationInvariant
            this = this@KerMorObject;
            this.registerProps('G');
        end
        
        function K = evaluate(this, x, y)
            % Evaluates the rotation and translation invariant kernel.
            %
            % Automatically yields the implementation of the
            % BaseKernel.evaluate function as only scalar evaluation is
            % required.
            sx = this.G*x;
            n1sq = sum(x.*sx,1);
            n1 = size(x,2);
            if nargin == 2;
                n2sq = n1sq;
                n2 = n1;
                y = x;
            else
                n2sq = sum(y.*(this.G*y),1);
                n2 = size(y,2);
            end;
            % Unfortunately here the square root has to be taken, as the
            % argument is simply ||a-b|| and not ||a-b||^2. the gaussian
            % overrides the evaluate fcn to avoid unnecessary rooting and
            % squaring.
            K = sqrt((ones(n2,1)*n1sq)' + ones(n1,1)*n2sq - 2*sx'*y);
            K = this.evaluateScalar(K);
        end
    end
    
    methods(Abstract)
        % Allows the evaluation of the function `\phi(x)` for scalar `x`
        % directly.
        Kx = evaluateScalar(x);
    end
    
end

