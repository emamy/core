classdef AKernelCoreFun < dscomponents.ACoreFun
    %KERNELCOREFUN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The function that combines the sub (time/system/param) kernels.
        % Must be a function handle that takes three arguments.
        SubKernelCombinationFun = @(t,s,p)t .* s .* p;
        
        % The Kernel to use for time variables
        %
        % See also: SystemKernel ParamKernel
        TimeKernel;
        
        % The Kernel to use for system variables
        %
        % See also: TimeKernel ParamKernel
        SystemKernel;
        
        % The Kernel to use for parameter variables
        %
        % See also: TimeKernel SystemKernel
        ParamKernel;
        
        % The state variable Snapshots used in the approximation.
        %
        % This is the union of all center data used within the kernel function.
        snData;
    end
    
    properties(SetAccess=private, Dependent)
        % A flag that tells subclasses whether the approximating kernel
        % function is rotation invariant. Depending on the flag projection
        % of subclasses can be different.
        RotationInvariantKernel;
    end
    
    methods
        function this = AKernelCoreFun
            % Default constructors. Initializes default kernels.
            this.TimeKernel = kernels.LinearKernel;
            this.SystemKernel = kernels.GaussKernel(2);
            % The default kernel is a neutral (=1) kernel as not all models
            % have parameters.
            this.ParamKernel = kernels.NoKernel;
            
            this.CustomProjection = true;
        end
        
        function target = project(this, V, W, target)
            if this.RotationInvariantKernel
                target.snData.xi = W' * this.snData.xi;
            end
        end
        
        function K = evaluateKernel(this, x, t, mu)
            % Evaluates the kernel using the given parameters both as first
            % and second argument.
            K = this.SubKernelCombinationFun(...
                this.TimeKernel.evaluate(t), ...
                this.SystemKernel.evaluate(x), ...
                this.ParamKernel.evaluate(mu));
        end
        
        function K = evaluateAtCenters(this, x, t, mu)
            % Evaluates the specified Kernel using the Sub-Kernel
            % combination function.
            %
            % Parameters:
            % x: The first set of state variables
            % t: The associated time steps of x (if used, take [] else)
            % mu: The parameters associated with x (if used, take []
            % else)
            %
            % @todo Ggf. Cache-Funktion wieder einbauen (feasibility?)
            % @todo check warning in this function.
            
            d = this.snData;
%             if (isempty(d.mui) && ~isa(this.ParamKernel,'kernels.NoKernel'))
%                 error('ad:ADG','mui is empty! ..');
%             end
            K = this.SubKernelCombinationFun(...
                this.TimeKernel.evaluate(d.ti, t), ...
                this.SystemKernel.evaluate(d.xi, x), ...
                this.ParamKernel.evaluate(d.mui, mu));
        end
        
        function value = get.RotationInvariantKernel(this)
            value = this.TimeKernel.RotationInvariant && ...
                this.SystemKernel.RotationInvariant && ...
                this.ParamKernel.RotationInvariant;
        end
    end
    
    methods(Access=protected)
        function target = clone(this, target)
            % Copy local variables
            target.snData = this.snData;
            target.SubKernelCombinationFun = this.SubKernelCombinationFun;
            
            % @todo Check whether kernels should be deepcopied, too
            target.TimeKernel = this.TimeKernel;
            target.SystemKernel = this.SystemKernel;
            target.ParamKernel = this.ParamKernel;
        end
    end
end

