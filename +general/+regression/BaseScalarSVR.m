classdef BaseScalarSVR < KerMorObject & ICloneable & IKernelCoeffComp
    %SCALARSVR Scalar support vector regression.
    %
    % Base class for any scalar SVR algorithm.
    %
    % @author Daniel Wirtz @date 2010-03-11
    %
    % @change{0,5,dw,2011-11-09} Also allowing to pass a double matrix to the setter for the K
    % property.
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
    % @todo
    % - Implement hyperparameter estimations from [CM04]!
    % - Look into the "pattern search" algorithm from Momma&Bennet 2002
    % - Create single MaxIterations property for all SVRs (if applicable)
    
    properties(SetObservable)
        % Minimum relative value for any alpha to be considered a support vector
        % coefficient.
        %
        % Relative means in this context alpha values divided by the absolute maximum over all values.
        %
        % @propclass{optional} The threshold of 16 magnitudes below the maximum coefficient should
        % be small enough to not have any bad influence but improve performance.
        %
        % @default: eps (2.2e-16) @type double
        AlphaRelMinValue = eps;
    end
    
    properties(Transient)
        % The kernel matrix to use.
        %
        % The reason why this is a property and not an argument is that
        % once a matrix is set multiple regressions for the same base
        % vector set can be performed easily.
        %
        % @propclass{data} Needed for the SVR to run in the first place.
        %
        % @type data.FileMatrix
        K;
    end
    
    properties(SetObservable, Dependent)
        % The regularization parameter `\lambda` for the primary minimization problem.
        %
        % @propclass{critical} Overly regularized functions may not approximate the data correctly,
        % while small `\lambda` lead to high coefficient values.
        %
        % @default 1 @type double
        Lambda;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        % The weighting of the slack variables.
        % 
        % Gets computed when Lambda is set, equals `C =
        % \frac{1}{2\lambda}`.
        %
        % @type double
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
            if (isa(value,'handle') && ~isa(value, 'data.FileMatrix')) || (~ismatrix(value) || ~isa(value,'double'))
                error('Value must be a data.FileMatrix or a double matrix.');
            end
            this.K = value;
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
        
        %% IKernelCoeffComp interface members
        function init(this, kexp)
            % Sets the kernel matrix.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            this.K = kexp.getKernelMatrix;
        end
        
        function [ci, svidx, sf] = computeKernelCoefficients(this, yi, initialai)
            % Implementation of the kernels.ICoeffComp interface
            %
            % Parameters:
            % yi: The target values `y_i` as row vector @type rowvec
            % initialai: Initial values for the coefficients `c_i`
            %
            % Return values:
            % ci: The coefficients `c_i` as row vector @type rowvec
            % svidx: The support vector indices of all elements of `c_i`
            % that regarded to be support vectors. @type integer
            %
            % See also: AlphaRelMinValue
            [ci, sf] = this.regress(yi, initialai);
            svidx = find(abs(ci) ./ max(abs(ci)) > this.AlphaRelMinValue);
            ci = ci(svidx);
            if isempty(svidx)
                ci = 0;
                svidx = 1;
                sf = StopFlag.NO_SUPPORT_VECTORS_FOUND;
            end
        end
    end
        
    methods(Abstract)
        % Performs the actual regression (template method)
        %
        % Parameters:
        % fxi: The `f(x_i)` values to regress given the kernel matrix `K` computed from
        % `\Phi(x_i,x_j)`.
        % initialai: The initial values to use for the coefficients `c_i`.
        % @default [] @type rowvec
        %
        % Return values:
        % ci: The kernel expansion coefficients `c_i`.
        % sf: The stop flag @type integer
        %
        % See also: StopFlag
        [ci, sf] = regress(this, fxi, initialai);
    end
    
    methods(Static)
        function c = getDefaultConfig
            c = general.regression.EpsSVRConfig([1; 1]);
        end
    end
    
end