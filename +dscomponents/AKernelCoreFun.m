classdef AKernelCoreFun < dscomponents.ACoreFun
    % Base class for kernel based system core functions
    % @author Daniel Wirtz @date April 2010
    %
    % This class is the composition of three different kernels, which are
    % one for time, parameters and state variables. Those kernels are
    % combined using the function handle set by the property
    % SubKernelCombinationFun.
    %
    % @change{0,2,dw,2011-03-21}
    % - Moved the dscomponents.AKernelCoreFun.RotationInvariant property
    % into a local private variable that gets updated only when the kernels
    % are changed. This way a repetitive in-simulation evaluation of the
    % property gets avoided.
    % - Changed & documented the default kernels.
    % - Added set method for
    % dscomponents.AKernelCoreFun.SubKernelCombinationFun property.
    
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
        
        % The state variable data used in the approximation.
        %
        % This is the union of all center data used within the kernel
        % function as a struct with the following fields: 
        % xi: The state variable centers
        % ti: The times at which the state xi is taken
        % mui: The parameter used to obtain the state xi
        snData;
    end
    
    properties(SetAccess=private)
        % A flag that tells subclasses if this (kernel-based)
        % function is rotation invariant.
        % 
        % Depending on the flag projection of subclasses can be different.
        RotationInvariant = false;
    end
    
    methods
        function this = AKernelCoreFun
            % Default constructors. Initializes default kernels.
            this.TimeKernel = kernels.NoKernel;
            this.SystemKernel = kernels.GaussKernel(10);
            % The default kernel is a neutral (=1) kernel as not all models
            % have parameters.
            this.ParamKernel = kernels.NoKernel;            
        end
        
        function target = project(this, V, W, target)
            % Call superclass method
            target = project@dscomponents.ACoreFun(this,V,W,target);
            % For rotation invariant kernel expansions the snapshots can be
            % transferred into the subspace without loss.
            if this.RotationInvariant
                target.snData.xi = W' * this.snData.xi;
            end
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
                        
            V = 1;
            if ~this.RotationInvariant && ~isempty(this.V)
                V = this.V;
            end
            K = this.SubKernelCombinationFun(...
                this.TimeKernel.evaluate(this.snData.ti, t), ...
                this.SystemKernel.evaluate(this.snData.xi, V*x), ...
                this.ParamKernel.evaluate(this.snData.mui, mu));
            % Checked alternative: not much faster so more general method
            % preferrable
%             K = this.TimeKernel.evaluate(this.snData.ti, t).* ...
%                 this.SystemKernel.evaluate(this.snData.xi, V*x).* ...
%                 this.ParamKernel.evaluate(this.snData.mui, mu);
        end
        
        function target = clone(this, target)
            % Copy local variables
            target.snData = this.snData;
            target.SubKernelCombinationFun = this.SubKernelCombinationFun;
            
            % @todo Check whether kernels should be deepcopied, too
            target.TimeKernel = this.TimeKernel;
            target.SystemKernel = this.SystemKernel;
            target.ParamKernel = this.ParamKernel;
        end
        
        function set.SubKernelCombinationFun(this, fhandle)
            if ~isa(fhandle,'function_handle')
                error('SubKernelCombinationFun must be a function handle.');
            end
            this.SubKernelCombinationFun = fhandle;
        end
    
        
        function set.ParamKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.ParamKernel = value;
                this.updateRotInv;%#ok
            else
                error('ParamKernel must be a subclass of kernels.BaseKernel.');
            end
        end
        
        function set.SystemKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.SystemKernel = value;
                this.updateRotInv;%#ok
            else
                error('SystemKernel must be a subclass of kernels.BaseKernel.');
            end
        end
        
        function set.TimeKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.TimeKernel = value;
                this.updateRotInv;%#ok
            else
                error('TimeKernel must be a subclass of kernels.BaseKernel.');
            end
        end
    end
    
    methods(Access=private)
        function updateRotInv(this)
            % Updates the RotationInvariant property of this CoreFun by
            % checking all registered kernels.
            this.RotationInvariant = true;
            if ~isempty(this.TimeKernel)
                this.RotationInvariant = this.RotationInvariant ...
                    && this.TimeKernel.RotationInvariant;
            end
            if ~isempty(this.SystemKernel)
                this.RotationInvariant = this.RotationInvariant ...
                    && this.SystemKernel.RotationInvariant;
            end
            if ~isempty(this.ParamKernel)
                this.RotationInvariant = this.RotationInvariant ...
                    && this.ParamKernel.RotationInvariant;
            end
        end
    end
end

