classdef BaseScalarSVR < KerMorObject & ICloneable & approx.IKernelCoeffComp
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
        %
        % @propclass{data} Needed for the SVR to run in the first place.
        K;    
    end
    
    properties(SetObservable)
        % The weighting of the slack variables. For larger C, the slack
        % variables are forced towards zero so that violations of the
        % eps-tube are getting penalized harder.
        % Theoretically, for C=Inf all fxi must be inside the eps-tube
        % around the original function.
        %
        % @propclass{critical} Large `C` leads to possibly large coefficients. Strongly linked to
        % the `nu` property of the ScalarNuSVR and the `\epsilon` property of the ScalarEpsSVR.
        %
        % @default 5
        %
        % See also: ScalarNuSVR ScalarEpsSVR
        C = 5;
        
        % Minimum relative value for any alpha to be considered a support vector
        % coefficient.
        %
        % Relative means in this context alpha values divided by the absolute maximum over all values.
        %
        % @propclass{optional} The threshold of 16 magnitudes below the maximum coefficient should
        % be small enough to not have any bad influence but improve performance.
        %
        % @default: eps (2.2e-16)
        AlphaRelMinValue = eps;
        
        % The quadratic solver internally used
        %
        % @propclass{optional} Different solvers should have different performance but should not
        % change the result. qpOASES so far is the fastest solver available in KerMor.
        %
        % @default solvers.qp.qpOASES
        % @type solvers.qp.BaseQPSolver
        QPSolver;
    end
    
    methods
        function this = BaseScalarSVR
            this = this@KerMorObject;
            this.QPSolver = solvers.qp.qpOASES;
            
            this.registerProps('K','C','AlphaRelMinValue','QPSolver');
        end
        
        function target = clone(this, target)
            target.K = this.K;
            target.C = this.C;
            target.AlphaRelMinValue = this.AlphaRelMinValue;
            target.QPSolver = this.QPSolver.clone;
        end
        
        function set.K(this, value)
            % Make matrix symmetric (can be false due to rounding errors)
            this.K = .5*(value + value');
        end
        
        %% approx.IKernelCoeffComp interface members
        function init(this, K)
            this.K = K;
        end
        
        function [ai, svidx] = computeKernelCoefficients(this, yi)
            [ai, svidx] = this.regress(yi);
            if isempty(svidx) 
                if any(yi ~= 0)
                    error('No support vectors found. Problem unsolvable with current config?\nQuadprog exit flag:%d\nQuadprog out.message:%s',exitflag,out.message);
                else
                    % Otherwise "fake" a support vector with zero coefficient!
                    ai = 0;
                    svidx = 1;
                end
            end
        end
    end
        
    methods(Abstract)
        [ai, svidx] = regress(this, fxi);
    end
    
end