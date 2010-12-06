classdef BaseScalarSVR < ICloneable
    %SCALARSVR Scalar support vector regression.
    %
    % Base class for any scalar SVR algorithm.
    %
    % @todo: Implement hyperparameter estimations from [CM04]!
    % @todo: Look into the "pattern search" algorithm from Momma&Bennet
    % 2002
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    properties
        % The kernel matrix to use.
        %
        % The reason why this is a property and not an argument is that
        % once a matrix is set multiple regressions for the same base
        % vector set can be performed easily.
        K;
        
        % The weighting of the slack variables. For larger C, the slack
        % variables are forced towards zero so that violations of the
        % eps-tube are getting penalized harder.
        % Theoretically, for C=Inf all fxi must be inside the eps-tube
        % around the original function.
        %
        % See also: eps
        C = 10;
        
        % Minimum value for any gamma to be considered a support vector
        % coefficient
        %
        % Default: 10*eps
        AlphaMinValue = 10*eps;
        
        % The quadratic solver internally used
        %
        % Defaults to solvers.qpOASES
        % @type solvers.BaseQPSolver
        QPSolver;
    end
    
    methods
        function this = BaseScalarSVR
            this.QPSolver = solvers.qpOASES;
        end
        
        function target = clone(this, target)
            target.K = this.K;
            target.C = this.C;
            target.AlphaMinValue = this.AlphaMinValue;
            target.QPSolver = this.QPSolver.clone;
        end
        
        function set.K(this, value)
            % Make matrix symmetric (can be false due to rounding errors)
            this.K = .5*(value + value');
        end
    end
        
    methods(Abstract)
        [ai,b,svidx] = regress(this, fxi);
    end
    
end

