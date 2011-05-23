classdef CompwiseKernelCoreFun < dscomponents.AKernelCoreFun & ...
        dscomponents.IGlobalLipschitz & dscomponents.IJacobian
    % Base class for component wise kernel expansion core functions.
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
    % @todo IMPORTANT: Allow different kernels for different centers!
    
    properties
        % The coefficient data for each dimension.
        %
        % @propclass{data}
        Ma;
    end
    
    properties(SetAccess=private, Dependent)
        Ma_norms;
    end
    
    methods
        
        function this = CompwiseKernelCoreFun
            this = this@dscomponents.AKernelCoreFun;
            
            this.CustomProjection = true;
            
            this.registerProps('Ma');
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            fx = this.Ma * this.evaluateAtCenters(x, t, mu)';
        end
        
        function J = getStateJacobian(this, x, t, mu)
            % Evaluates the jacobian matrix of this function at the given point.
            %
            % As this is a component-wise kernel expansion, the jacobian is easily computed using
            % the coefficient vectors and the state variable nablas.
            %
            % See also: dscomponents.IJacobian kernels.BaseKernel
            N = this.evaluateStateNabla(x, t, mu);
            J = this.Ma * N';
        end
        
        function c = getGlobalLipschitz(this, t, mu)%#ok
            % @todo validate computation
            %abs(this.TimeKernel.evaluate(this.Centers.ti,t).*this.ParamKernel.evaluate(this.Centers.mui,mu));
            c = sum(this.Ma_norms) * this.SystemKernel.getGlobalLipschitz;
            warning('KerMor:globallipschitz','not yet implemented/validated correctly!');
        end
        
        function target = clone(this, target)
            % Check if already a subclass instance is passed upwards and
            % use it if given; otherwise create a new instance here!
            if nargin < 2
                target = dscomponents.CompwiseKernelCoreFun;
            end
            
            target = clone@dscomponents.AKernelCoreFun(this, target);
            
            % Copy local variables
            target.Ma = this.Ma;
        end
        
        function m = get.Ma_norms(this)
            m = sqrt(sum(this.Ma.^2,1));
        end
    end
        
    methods(Sealed)
        
        function copy = project(this, V, W)
            % IMPORTANT NOTE:
            % This function was put into the sealed methods block
            % since the project function does not accept the 4th argument
            % "target". This is due to the fact that no further subclasses
            % are intended for this class; if this gets the case, one needs
            % to add the target argument and check if it's set or not:
            % @code
            %if nargin < 4
            %    copy = this.clone;
            %else
            %    copy = target;
            %end
            % @endcode
            
            % Create copy
            copy = this.clone;
            
            % Call superclass project with new instance
            copy = project@dscomponents.AKernelCoreFun(this, V, W, copy);
            
            % Project coefficients & offsets
            copy.Ma = W' * this.Ma;
        end
    end
    
end

