classdef ParamTimeKernelExpansion < kernels.KernelExpansion
% ParamTimeKernelExpansion: Kernel expansion class for time and/or parameter dependent kernels.
%
% This class represents a function `f` in the space induced by given kernels `\Phi_{s,t,\mu}` for state, time
% and parameter space as ``f(x,t,\mu) = \sum\limits_{i=1}^N c_i\Phi_s(x,x_i)\Phi_t(t,t_i)\Phi_\mu(\mu,\mu_i)``
%
% This class is the composition of three different kernels, which are
% one for time, parameters and state variables. Those kernels are
% combined using the function handle set by the property
% SubKernelCombinationFun.
%
% @author Daniel Wirtz @date 2011-07-07
%
% @new{0,5,dw,2011-07-07} Added this class.
%
% Former changes from old class AKernelCoreFun:
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @change{0,3,sa,2011-04-16} Implemented Setter for the property 'Centers'
%
% @new{0,3,dw,2011-04-15} Added the dscomponents.AKernelCoreFun.evaluateStateNabla method to
% allow efficient computation of kernel expansion jacobians.
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
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable, Dependent)
        % The Kernel to use for time variables
        %
        % @propclass{critical} Correct choice of the time kernel greatly influences the function
        % behaviour.
        %
        % @default kernels.NoKernel
        %
        % See also: Kernel ParamKernel
        TimeKernel;
        
        % The Kernel to use for parameter variables
        %
        % @propclass{critical} Correct choice of the system kernel greatly influences the function
        % behaviour.
        %
        % @default kernels.NoKernel
        %
        % See also: TimeKernel Kernel
        ParamKernel;
    end
    
    properties(SetObservable)
        % The function that combines the sub (time/system/param) kernels.
        % Must be a function handle that takes three arguments.
        %
        % @propclass{optional} Standard kernel combination is pointwise multiplication.
        %
        % @default @code @(t,s,p)t .* s .* p @endcode
        %
        % See also: StateNablaCombinationFun
        SubKernelCombinationFun = @(t,s,p)t .* s .* p;
        
        % The combination function for the nabla of the system/state kernel with the other kernels.
        %
        % Depends directly on the SubKernelCombinationFun
        %
        % @propclass{optional} Adjust this whenever the SubKernelCombinationFun is changed.
        %
        % See also: SubKernelCombinationFun
        StateNablaCombinationFun = @(t,s,p)bsxfun(@times, s, t.*p);
    end
    
    properties(Access=private)
        fTK;
        fPK;
    end
    
    methods
        function this = ParamTimeKernelExpansion
            % Default constructor. 
            %
            % Initializes default kernels.
            
            % Set custom projection to true as the project method is
            % overridden
            this = this@kernels.KernelExpansion;
            

            
%             this.fSK = kernels.GaussKernel;
            % The default kernels for time and parameters are neutral (=1)
            % kernels as not all models have time or parameter dependent
            % system functions.
            this.fTK = kernels.NoKernel;
            this.fPK = kernels.NoKernel;
            
            % DONT CHANGE THIS LINE unless you know what you are doing or
            % you are me :-)
            this.updateRotInv;
            
            this.registerProps('TimeKernel','ParamKernel',...
                'SubKernelCombinationFun','StateNablaCombinationFun');
        end

        function fx = evaluate(this, x, t, mu)
%             phi = this.SubKernelCombinationFun(...
%                 this.fTK.evaluate(t, this.Centers.ti), ...
%                 this.fSK.evaluate(x, this.Centers.xi), ...
%                 this.fPK.evaluate(mu, this.Centers.mui));
            fx = this.Ma * this.getKernelVector(x, t, mu)';
        end
        
        function phi = getKernelVector(this, x, t, mu)
            % Returns the kernel vector `\varphi(x,t,\mu) = (\Phi(x,x_i)\Phi(t,t_i)\Phi(\mu,\mu_i))_i`.
            % (for the case of a product 'SubKernelCombinationFun'
            phi = this.SubKernelCombinationFun(...
                this.fTK.evaluate(t, this.Centers.ti), ...
                this.fSK.evaluate(x, this.Centers.xi), ...
                this.fPK.evaluate(mu, this.Centers.mui));
        end
                
        function J = getStateJacobian(this, x)
            % Evaluates the jacobian matrix of this function at the given point.
            %
            % As this is a component-wise kernel expansion, the jacobian is easily computed using
            % the coefficient vectors and the state variable nablas.
            %
            % See also: dscomponents.IJacobian kernels.BaseKernel
            if size(x, 2) > 1
                error('Derivaties only possible for single vector/point, as already returning a matrix with derivatives at all centers.');
            end
            N = this.StateNablaCombinationFun(...
                this.fTK.evaluate(t, this.Centers.ti),...
                this.fSK.getNabla(x, this.Centers.xi),...
                this.fPK.evaluate(mu, this.Centers.mui));
            J = this.Ma * N';
        end
        
        function c = getGlobalLipschitz(this, t, mu)
            % Overrides the implementation in KernelExpansion.
            % @todo validate computation
            k = this.TimeKernel.evaluate(this.Centers.ti,t).*this.ParamKernel.evaluate(this.Centers.mui,mu);
            c = sum(this.Ma_norms .* k) * this.Kernel.getGlobalLipschitz;
            %warning('KerMor:globallipschitz','not yet implemented/validated correctly!');
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
                getKernelMatrix@kernels.KernelExpansion(this), ...
                this.fPK.evaluate(this.Centers.mui));
        end
        
        function copy = clone(this, copy)
            if nargin == 1
                copy = kernels.ParamTimeKernelExpansion;
            end
            % Copy local variables
            copy = clone@kernels.KernelExpansion(this, copy);
            copy.SubKernelCombinationFun = this.SubKernelCombinationFun;
            copy.StateNablaCombinationFun = this.StateNablaCombinationFun;
            copy.fTK = this.fTK;
            copy.fPK = this.fPK;
        end       
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
        
        function set.TimeKernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.fTK = value;
                this.updateRotInv;
            else
                error('TimeKernel must be a subclass of kernels.BaseKernel.');
            end
        end
    end
    
    methods(Access=protected)
        function updateRotInv(this)
            % Updates the RotationInvariant property of this KernelExpansion
            updateRotInv@kernels.KernelExpansion(this);
            this.RotationInvariant = this.RotationInvariant && ...
                isa(this.fTK,'kernels.IRotationInvariant') && ...
                isa(this.fPK,'kernels.IRotationInvariant');
        end
    end
    
end