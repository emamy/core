classdef BaseKernelApprox < approx.BaseApprox
    %BASEKERNELAPPROX Summary of this class goes here
    %   Detailed explanation goes here
    % IMPORTANT:
    % Changing the Kernels AFTER a model's offline simulation but BEFORE
    % the projection process (aka buildReducedModel)
    % a way that the rotation invariance changes may lead to inefficient projection
    
    properties
        % The function that combines the sub (time/system/param) kernels.
        % Must be a function handle that takes three arguments.
        SubKernelCombinationFun = @(t,s,p)t .* s .* p;
        
        % The Kernel to use for time variables
        %
        % See also: SystemKernel ParamKernel
        TimeKernel = kernels.RBFKernel;
        
        % The Kernel to use for system variables
        %
        % See also: TimeKernel ParamKernel
        SystemKernel = kernels.RBFKernel;
        
        % The Kernel to use for parameter variables
        %
        % See also: TimeKernel SystemKernel
        ParamKernel = kernels.RBFKernel;
    end
    
    properties(SetAccess=private, GetAccess=protected, Dependent)
        % A flag that tells subclasses whether the approximating kernel
        % function is rotation invariant. Depending on the flag projection
        % of subclasses can be different.
        RotationInvariantKernel;
    end
    
    methods(Access=protected)
        function K = evaluateKernel(this, x, y)
            % Evaluates the specified Kernel by using a multiplication
            % kernel (default with weights=1 until application leads to it)
            %
            % Argument y is optional, if not given y=x is assumed.
            
            [t1,x1,mu1] = this.splitTripleVect(x);
            
            % NOTE:
            % The "weird" structure of this code is for speed reasons.
            % The computation of the kernel matrix will be of a greater
            % magnitude than some if-clauses; important (can) be the
            % correct call of the kernels.evaluate methods as they may be
            % more efficient for x=y
            
            if nargin == 3
                [t2,x2,mu2] = this.splitTripleVect(y);
            end
            
            % Make sure the ParamKernel doesnt destroy the other kernels in
            % case no parameters are given. so set to 1 per default. Also
            % for speed, dont bother to call evaluate on param kernels if
            % empty.
            pker = 1;
            if ~isempty(mu1)
                if nargin == 2
                    pker = this.ParamKernel.evaluate(mu1);
                else
                    pker = this.ParamKernel.evaluate(mu1,mu2);
                end
            end
            % For speed reasons call the reduced version of the kernel if
            % y=x holds
            if nargin == 2
                K = this.SubKernelCombinationFun(...
                    this.TimeKernel.evaluate(t1), ...
                    this.SystemKernel.evaluate(x1), ...
                    pker);
            else
                K = this.SubKernelCombinationFun(...
                    this.TimeKernel.evaluate(t1,t2), ...
                    this.SystemKernel.evaluate(x1,x2), ...
                    pker);
            end
        end
        
        function target = clone(this, target)
            % Clones this instance into the specified target (must be
            % subclass) or creates a new instance.
            if nargin == 1
                target = BaseKernelApprox;
            elseif ~isa(target,'approx.BaseKernelApprox')
                error('Cloning into a non-approx.BaseKernelApprox instance impossible.');
            end
            % Call superclass' clone method
            clone@approx.BaseApprox(this, target);
            % Copy local values into new object
            target.SubKernelCombinationFun = this.SubKernelCombinationFun;
            target.TimeKernel = this.TimeKernel;
            target.SystemKernel = this.SystemKernel;
            target.ParamKernel = this.ParamKernel;
            % TODO: Check whether kernels should be deepcopied, too
        end
    end
    
    methods
        function value = get.RotationInvariantKernel(this)
            value = this.TimeKernel.RotationInvariant && ...
                this.SystemKernel.RotationInvariant && ...
                this.ParamKernel.RotationInvariant;
        end
        
        % @TODO: getter & setter
    end
    
end

