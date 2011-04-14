classdef AKernelCoreFun < dscomponents.ACoreFun
    % Base class for kernel based system core functions
    % @author Daniel Wirtz @date April 2010
    %
    % This class is the composition of three different kernels, which are
    % one for time, parameters and state variables. Those kernels are
    % combined using the function handle set by the property
    % SubKernelCombinationFun.
    %
    % @change{0,3,dw,2011-04-6} Set the kernel properties to dependent and
    % introduced private storage members for them. Also improved the setter
    % for the SubKernelCombinationFun (checks for 3 args now)
    %
    % @new{0,3,dw,2011-04-04} Added the
    % dscomponents.AKernelCoreFun.getKernelMatrix method
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
        
        % The kernel centers used in the approximation.
        %
        % This is the union of all center data used within the kernel
        % function as a struct 
        % with the following fields:
        % - 'xi' The state variable centers
        % - 'ti' The times at which the state xi is taken
        % - 'mui' The parameter used to obtain the state xi
        %
        % The only required field is xi, others can be set to [] if not
        % used.
        %
        % @default @code struct('xi',[],'ti',[],'mui',[]) @endcode
        Centers;
    end
    
    properties(Dependent)
        % The Kernel to use for time variables
        %
        % @default kernels.NoKernel
        % See also: SystemKernel ParamKernel
        TimeKernel;
        
        % The Kernel to use for system variables
        %
        % @default kernels.GaussKernel
        % See also: TimeKernel ParamKernel
        SystemKernel;
        
        % The Kernel to use for parameter variables
        %
        % @default kernels.NoKernel
        % See also: TimeKernel SystemKernel
        ParamKernel;
    end
    
    properties(SetAccess=private)
        % A flag that tells subclasses if this (kernel-based)
        % function is rotation invariant.
        % 
        % Depending on the flag projection of subclasses can be different.
        RotationInvariant = false;
    end
    
    properties(Access=private)
        fSK;
        fTK;
        fPK;
    end
    
    methods
        function this = AKernelCoreFun
            % Default constructor. 
            %
            % Initializes default kernels.
            
            % Set custom projection to true as the project method is
            % overridden
            this.CustomProjection = true;
            % Kernel based core functions allow for multi-argument evaluations by nature.
            this.MultiArgumentEvaluations = true;
            
            this.fSK = kernels.GaussKernel(10);
            % The default kernels for time and parameters are neutral (=1)
            % kernels as not all models have time or parameter dependent
            % system functions.
            this.fTK = kernels.NoKernel;
            this.fPK = kernels.NoKernel;
            
            % DONT CHANGE THIS LINE unless you know what you are doing or
            % you are me :-)
            this.updateRotInv;
        end
        
        function target = project(this, V, W, target)
            % Call superclass method
            target = project@dscomponents.ACoreFun(this,V,W,target);
            % For rotation invariant kernel expansions the snapshots can be
            % transferred into the subspace without loss.
            if this.RotationInvariant
                target.Centers.xi = W' * this.Centers.xi;
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
            % See also: getKernelMatrix
            %
            % @todo Ggf. Cache-Funktion wieder einbauen (feasibility?)
                        
            V = 1;
            if ~this.RotationInvariant && ~isempty(this.V)
                V = this.V;
            end
            K = this.SubKernelCombinationFun(...
                this.fTK.evaluate(this.Centers.ti, t), ...
                this.fSK.evaluate(this.Centers.xi, V*x), ...
                this.fPK.evaluate(this.Centers.mui, mu));
            % Checked alternative: not much faster so more general method
            % preferrable
%             K = this.TimeKernel.evaluate(this.Centers.ti, t).* ...
%                 this.SystemKernel.evaluate(this.Centers.xi, V*x).* ...
%                 this.ParamKernel.evaluate(this.Centers.mui, mu);
        end
        
        function K = getKernelMatrix(this)
            % Computes the kernel matrix for the currently set center data.
            %
            % Convenience method, equal to call evaluateAtCenters with
            % passing the center data as arguments. However, this method is
            % (possibly, depends on kernel implementation) faster as it
            % calls the kernels kernels.BaseKernel.evaluate -method with
            % one argument.
            %
            % See also: evaluateAtCenters
            K = this.SubKernelCombinationFun(...
                this.fTK.evaluate(this.Centers.ti), ...
                this.fSK.evaluate(this.Centers.xi), ...
                this.fPK.evaluate(this.Centers.mui));
        end
        
        function target = clone(this, target)
            % Copy local variables
            target.Centers = this.Centers;
            target.SubKernelCombinationFun = this.SubKernelCombinationFun;
            
            % @todo Implement cloning for kernels too
            target.fTK = this.fTK;
            target.fSK = this.fSK;
            target.fPK = this.fPK;
            target.RotationInvariant = this.RotationInvariant;
        end
        
%         function autoConfGamma(this, model)
%             % Implements the template method from IAutoConfigure.
%             %
%             % For this class autoconfiguration means detection of the
%             % "ideal" radius for gaussian kernels, if used. The strategy is
%             % to enforce that for the largest distance between any two
%             % considered centers the sum of both nearby kernel evaluations
%             % equals one, i.e.
%             % ``e^{-\frac{\left(\frac{d}{2}\right)^2}{\gamma}} =
%             % \frac{1}{2}``
%             % if `d` is the largest distance.
%             % 
%             % Parameters:
%             % model: The current model instance
%             %
%             % See also: IAutoConfigure
%             
%             % Settings.
%             %zero = 1e-4;
%             %trange = 3; % nonzero over trange times the dt-distance
%             %srange = 3; % nonzero over srange times the maximum distance within the training data
%             %prange = 2; % nonzero over prange times the param samples distance
%             data = model.Data.ApproxTrainData;
%             v = unique(round(data(1,:)));
%             
%             %% State kernel gamma
%             if isa(this.SystemKernel,'kernels.GaussKernel')
%                 xd = sqrt(sum(data(4:end,:).^2));
% 
%                 % Find samples for each parameter    
%                 maxdiff = zeros(1,length(v));
%                 for muidx = 1:length(v)
%                     sel = data(1,:) == v(muidx);
%                     tmp = xd(sel);
%                     maxdiff(muidx) = max(abs(tmp(1:end-1)-tmp(2:end)));
%                 end
%                 d = max(maxdiff);
%                 this.SystemKernel.setGammaForDistance(d/2,.5);
%             end
%             
%             %% Time kernel
%             if isa(this.TimeKernel,'kernels.GaussKernel')
%                 warning('Code:unchecked','Implementation not yet finished/ideal!');
%                 this.TimeKernel.setGammaForDistance(model.dt/2,.5);
%             end
%             
%             %% Param kernel
%             if isa(this.ParamKernel,'kernels.GaussKernel')
%                 params = model.Data.getParams(v);
%                 % Gives a matrix with parameters in each column
%                 mud = sqrt(sum(params.^2));
%                 % Create distance matrix
%                 dist = abs(repmat(mud,size(mud,2),1)-repmat(mud',1,size(mud,2)));
%                 this.ParamKernel.setGammaForDistance(max(dist(:))/2,.5);
%             end
%         end
        
    end
    
    %% Getter & Setter
    methods    
        function set.SubKernelCombinationFun(this, fhandle)
            if ~isa(fhandle,'function_handle')
                error('SubKernelCombinationFun must be a function handle.');
            elseif nargin(fhandle) ~= 3
                error('SubKernelCombinationFun must take exactly three input arguments.');
            end
            this.SubKernelCombinationFun = fhandle;
        end
        
        function k = get.ParamKernel(this)
            k = this.fPK;
        end
        
        function k = get.SystemKernel(this)
            k = this.fSK;
        end
        
        function k = get.TimeKernel(this)
            k = this.fTK;
        end
    
        function set.ParamKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.fPK = value;
                this.updateRotInv;
            else
                error('ParamKernel must be a subclass of kernels.BaseKernel.');
            end
        end
        
        function set.SystemKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.fSK = value;
                this.updateRotInv;
            else
                error('SystemKernel must be a subclass of kernels.BaseKernel.');
            end
        end
        
        function set.TimeKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.fTK = value;
                this.updateRotInv;
            else
                error('TimeKernel must be a subclass of kernels.BaseKernel.');
            end
        end
    end
    
    methods(Access=private)
        function updateRotInv(this)
            % Updates the RotationInvariant property of this CoreFun by
            % checking all registered kernels.
            this.RotationInvariant = ...
                this.fTK.RotationInvariant &&...
                this.fSK.RotationInvariant && ...
                this.fPK.RotationInvariant;
        end
    end
end

