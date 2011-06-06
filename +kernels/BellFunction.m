classdef BellFunction < kernels.BaseKernel & kernels.IRotationInvariant
    %BELLFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @todo export different estimator strategies into external classes =>
    % simulation speedup, separation of concerns..
    %
    % @todo investigate why the newton iteration sometimes exceeds max
    % iteration limit (example: intermittently test_LinearModelParams is
    % such a case)
    %
    % @docupdate Properties and class description
    %
    % @change{0,4,dw,2011-05-27} Changed `x_0` to `r_0` and `x_R` to `r_m` as adopted in the WH10
    % Paper.
    %
    % @change{0,4,dw,2011-05-20} Removed the getLipschitzFunction method as it causes LARGE overhead
    % when being called very often. Instead, the estimator always uses the improved estimation
    % procedure.
    %
    % @change{0,4,dw,2011-05-19} 
    % - Removed the PenaltyFactor property as no longer needed.
    % - Modified the newton iteration penalization such that the
    % 2nd degree polynomials extend the objective function both in value and derivative. This
    % ensures propert continuation of the objective function into the penalized area with respect to
    % the current kernel configuration (i.e. kernels.GaussKernel: small Gamma property causes high
    % derivatives and thus large estimations.
    %
    % @new{0,3,dw,2011-04-06} Added a new property
    % kernels.BellFunction.MaxNewtonIterations that avoids computations to
    % hang if the newton iteration does not come to a hold. An error will
    % be thrown as finding the correct minima is necessary.
    
    properties(SetObservable)
        % Point of maximum first derivative on scalar evaluation.
        %
        % @propclass{critical} This value is essential for any bell function.
        r0;
    end
    
    properties(SetAccess=private, Dependent)
        % The maximum ("right") value for any `r_s`.
        rm;
    end
    
    properties(Access=private, Transient)
        priv_rr = [];
    end
    
    methods
        
        function this = BellFunction
            this = this@kernels.BaseKernel;
            this.registerProps('r0');
        end
        
        function c = getGlobalLipschitz(this)
            % Computes the absolute value of the first derivative at x0
            % Implements the template method from BaseKernel.
            c = abs(this.evaluateD1(this.r0));
        end
        
        function set.r0(this, value)
            if ~isposrealscalar(value) 
                error('r0 must be a scalar greater than zero.');
            end
            this.r0 = value;
            this.priv_rr = [];%#ok
        end
                
        function value = get.rm(this)
            if isempty(this.priv_rr)
                this.priv_rr = this.evaluateScalar(0)*this.r0 /...
                    (this.evaluateScalar(0)-this.evaluateScalar(this.r0));
            end
            value = this.priv_rr;
        end
    end
    
    methods(Abstract)
        % Method for first derivative evaluation
        dr = evaluateD1(r);
        
        % Method for second derivative evaluation
        ddr = evaluateD2(r);
    end
    
end

