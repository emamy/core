classdef BaseScalarSVR < KerMorObject & ICloneable & approx.algorithms.IKernelCoeffComp
    %SCALARSVR Scalar support vector regression.
    %
    % Base class for any scalar SVR algorithm.
    %
    % @author Daniel Wirtz @date 11.03.2010
    %
    % @change{0,5,dw,2011-08-22} Added the regularization parameter Lambda and made the C constraint
    % dependent on that, as also done in literature. Moved the QPSolver to
    % the classes that really use one.
    %
    % @change{0,4,dw,2011-05-31} Removed the 'svidx' parameter from the regress interface.
    %
    % @change{0,3,sa,2011-05-07} Implemented Setter for the properties of
    % this class
    %
    % @todo: Implement hyperparameter estimations from [CM04]!
    % @todo: Look into the "pattern search" algorithm from Momma&Bennet
    % 2002
    
    properties(SetObservable)
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
    end
    
    properties
        % The kernel matrix to use.
        %
        % The reason why this is a property and not an argument is that
        % once a matrix is set multiple regressions for the same base
        % vector set can be performed easily.
        %
        % @propclass{data} Needed for the SVR to run in the first place.
        %
        % @type data.AKernelMatrix
        K;        
    end
    
    properties(SetObservable, Dependent)
        % The regularization parameter `\lambda` for the primary minimization problem.
        %
        % @propclass{critical} Overly regularized functions may not approximate the data correctly,
        % while small `\lambda` lead to high coefficient values.
        %
        % @default 1
        Lambda;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        % The weighting of the slack variables.
        % 
        % Gets computed when Lambda is set, equals `C =
        % \frac{1}{2\lambda}`.
        %
        % See also: ScalarNuSVR ScalarEpsSVR
        C = .5;
    end
        
    properties(Access=private)
        % The internal real value for lambda.
        fLambda = 1;
    end
    
    methods
        function this = BaseScalarSVR
            this = this@KerMorObject;
            this.registerProps('K','Lambda','AlphaRelMinValue');
        end
        
        function target = clone(this, target)
            target.K = this.K;
            target.C = this.C;
            target.AlphaRelMinValue = this.AlphaRelMinValue;
        end
        
        function set.K(this, value)
            if ~isa(value, 'data.AKernelMatrix')
                error('Value must be a data.AKernelMatrix');
            end
            % Make matrix symmetric (can be false due to rounding errors)
            this.K = value;
            %this.K = .5*(value + value');
        end
        
        function set.Lambda(this, value)
            if ~isposrealscalar(value)
                error('Lambda must be a positive real scalar');
            end
            this.fLambda = value;
            this.C = 1/(2*value);
        end
        
        function value = get.Lambda(this)
            value = this.fLambda;
        end
        
        function set.AlphaRelMinValue(this, value)
            if ~isposrealscalar(value)
                error('AlphaRelMinValue must be a positive real scalar');
            end
            this.AlphaRelMinValue = value;
        end
        
        %% approx.algorithms.IKernelCoeffComp interface members
        function init(this, K)
            this.K = K;
        end
        
        function [ai, svidx] = computeKernelCoefficients(this, yi, initialai)
            % @throws KerMor:coeffcomp:failed Thrown if the coefficient computation failed for a
            % reason known within KerMor.
            %
            % @throws KerMor:svr:nosupportvectors No support vectors could be found for the given
            % configuration.
            try
                ai = this.regress(yi, initialai);
                svidx = find(abs(ai) ./ max(abs(ai)) > this.AlphaRelMinValue);
                ai = ai(svidx);
            catch ME
                if strcmp(ME.identifier,'KerMor:solvers:qp:notconverged')
                    m = MException('KerMor:coeffcomp:failed','Regression failed');
                    m.addCause(ME);
                    m.throw;
                else
                    rethrow(ME);
                end
            end
            if isempty(svidx)
%                 if any(yi ~= 0)
%                     m = MException('KerMor:svr:nosupportvectors','No support vectors found. Problem unsolvable with current config?');
%                     m.throw;
%                 else
                    % Otherwise "fake" a support vector with zero coefficient!
                    ai = 0;
                    svidx = 1;
%                 end
            end
        end
    end
        
    methods(Abstract)
        % Performs the actual regression (template method)
        %
        % Parameters:
        % fxi: The `f(x_i)` values to regress given the kernel matrix `K` computed from
        % `\Phi(x_i,x_j)`.
        %
        % Return values:
        % ai: The kernel expansion coefficients `\alpha_i`.
        ai = regress(this, fxi, initialai);
    end
    
end