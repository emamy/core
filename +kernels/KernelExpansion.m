classdef KernelExpansion < KerMorObject & ICloneable & ...
        dscomponents.IJacobian & dscomponents.IGlobalLipschitz
% KernelExpansion: Base class for state-space kernel expansions
%
% The KernelExpansion class represents a function `f` in the space induced by a given kernel `\Phi`
% as ``f(x) = \sum\limits_{i=1}^N c_i\Phi(x,x_i)``
% 
% This class has been introduced in order to allow standalone use of kernel expansions in
% approximation contexts outside of the intended model reduction scheme.
%
% @author Daniel Wirtz @date 2011-07-07
%
% @change{0,5,dw,2011-07-28} Fixed the evaluate method, it had an argument Centers.xi too much.
%
% @new{0,5,dw,2011-07-07} Added this class.
%
% Changes from former class CompwiseKernelCoreFun:
%
% @change{0,4,dw,2011-05-03}
% - Removed the off property as proper RKHS functions dont have an offset
% - No more setter for Ma property - too slow even during offline phase
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @change{0,3,sa,2011-04-16} Implemented Setter for the property 'off'
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable, Dependent)
        % The Kernel to use for system variables
        %
        % @propclass{critical} Correct choice of the system kernel greatly influences the function
        % behaviour.
        %
        % @default kernels.GaussKernel
        % @type kernels.BaseKernel
        %
        % See also: kernels
        Kernel;
    end
    
    properties(SetAccess=private, Dependent)
        % The norms of the coefficient matrix of the kernel expansion
        %
        % @type rowvec
        Ma_norms;
    end
    
    properties(SetObservable)
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
        % @propclass{data}
        %
        % @default @code struct('xi',[],'ti',[],'mui',[]) @endcode
        % @type struct
        Centers;
        
        % The coefficient data for each dimension.
        %
        % @propclass{data}
        %
        % @type matrix
        Ma;
    end
    
    properties(SetAccess=protected)
        % A flag that tells subclasses if this (kernel-based)
        % function is rotation invariant.
        % 
        % Depending on the flag projection of subclasses can be different.
        %
        % @default false @type logical
        RotationInvariant = false;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        % The inner (state) kernel object
        %
        % @type kernels.BaseKernel @default kernels.GaussKernel;
        fSK;
    end
    
    methods
        function this = KernelExpansion
            % Default constructor. 
            %
            % Initializes default kernels.
            
            % Set custom projection to true as the project method is
            % overridden
            this = this@KerMorObject;
            
            this.fSK = kernels.GaussKernel;
            
            % DONT CHANGE THIS LINE unless you know what you are doing or
            % you are me :-)
            this.updateRotInv;
            
            this.registerProps('Kernel','Centers');
        end
        
        function fx = evaluate(this, x)
            %fx = this.Ma * this.fSK.evaluate(x, this.Centers.xi)';
            fx = this.Ma * this.getKernelVector(x)';
        end
        
        function phi = getKernelVector(this, x)
            % Returns the kernel vector `\varphi(x) = (\Phi(x,x_i))_i`.
            phi = this.fSK.evaluate(x, this.Centers.xi);
        end
                
        function J = getStateJacobian(this, x, varargin)
            % Evaluates the jacobian matrix of this function at the given point.
            %
            % As this is a component-wise kernel expansion, the jacobian is easily computed using
            % the coefficient vectors and the state variable nablas.
            %
            % Parameters:
            % x: The state vector at which the jacobian is required.
            % varargin: Dummy parameter to satisfy the interface calls for
            % kernels.ParamTimeKernelExpansion classes, which also give a
            % `t` and `\mu` parameter.
            %
            % See also: dscomponents.IJacobian kernels.BaseKernel
            if size(x, 2) > 1
                error('Derivaties only possible for single vector/point, as already returning a matrix with derivatives at all centers.');
            end
            J = this.Ma * this.fSK.getNabla(x, this.Centers.xi)';
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
            K = this.fSK.evaluate(this.Centers.xi);
        end
        
        function c = getGlobalLipschitz(this, t, mu)%#ok            
            c = sum(this.Ma_norms) * this.Kernel.getGlobalLipschitz;
        end
        
        function copy = clone(this, copy)
            if nargin == 1
                copy = kernels.KernelExpansion;
            end
            % Copy local variables
            copy.Centers = this.Centers;
            copy.Ma = this.Ma;
            copy.fSK = this.fSK;
            copy.RotationInvariant = this.RotationInvariant;
        end
        
        function m = get.Ma_norms(this)
            m = sqrt(sum(this.Ma.^2,1));
        end
    end
    
    %% Getter & Setter
    methods    
        function set.Centers(this, value)
            C = {'xi','ti','mui'};            
            if isempty(value) || any(isfield(value, C))
                this.Centers = value;
%                 if ~isfield(value,'xi')
%                     warning('REQUIRED_FIELD:Empty','xi is a required field, which is left empty');
%                 end
            else
                error('Value passed is not a valid struct. Only the field names xi, ti, mui are allowed');
            end            
        end
        
        function k = get.Kernel(this)
            k = this.fSK;
        end
    
        function set.Kernel(this, value)
            if isa(value,'kernels.BaseKernel')
                this.fSK = value;
                this.updateRotInv;
            else
                error('Kernel must be a subclass of kernels.BaseKernel.');
            end
        end
    end
    
    methods(Access=protected)
        function updateRotInv(this)
            % Updates the RotationInvariant property of this CoreFun by
            % checking all registered kernels.
            this.RotationInvariant = isa(this.fSK,'kernels.IRotationInvariant');
        end
    end
    
end