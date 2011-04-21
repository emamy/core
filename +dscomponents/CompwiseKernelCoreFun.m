classdef CompwiseKernelCoreFun < dscomponents.AKernelCoreFun & ...
        dscomponents.IGlobalLipschitz & dscomponents.IJacobian
    %COMPWISEKERNELCOREFUN Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @docupdate
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @todo IMPORTANT: Allow different kernels for different centers!
    %
    % @change{0,3,sa,2011-04-16} Implemented Setter for the property 'off'
    
    properties
        % The `b_k` offsets for each dimension.
        %
        % Set to empty if no offsets are used. This property is
        % preinitialized to [].
        %
        % @propclass{data}
        %
        % See also: Ma
        off = [];
    end
    
    properties(Dependent)
        % The coefficient data for each dimension.
        %
        % @propclass{data}
        Ma;
    end
    
    properties(Access=private)
        machanged = true;
        fMa;
        Ma_norms;
    end
    
    methods
        
        function this = CompwiseKernelCoreFun
            this = this@dscomponents.AKernelCoreFun;
            this.CustomProjection = true;
            
            this.registerProps('off','Ma');
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            % @todo: check for correct usage of the sets from kernel
            % expansions in error estimators! (not yet considered)
            if ~isempty(this.off)
                fx = this.Ma * this.evaluateAtCenters(x, t, mu)' + repmat(this.off,1,size(x,2));
            else
                fx = this.Ma * this.evaluateAtCenters(x, t, mu)';
            end
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
        
        function c = getGlobalLipschitz(this, t, mu)
            % @todo validate computation
            
            %k = abs(this.TimeKernel.evaluate(this.Centers.ti,t).*this.ParamKernel.evaluate(this.Centers.mui,mu));
            if this.machanged
                this.Ma_norms = sqrt(sum(this.Ma.^2,1));
                this.machanged = false;
            end
            c = sum(this.Ma_norms) * this.SystemKernel.getGlobalLipschitz;
            %warning('some:id','not yet implemented/validated correctly!');
        end
        
        function target = clone(this, target)
            % Check if already a subclass instance is passed upwards and
            % use it if given; otherwise create a new instance here!
            if nargin < 2
                target = dscomponents.CompwiseKernelCoreFun;
            end
            
            target = clone@dscomponents.AKernelCoreFun(this, target);
            
            % Copy local variables
            target.fMa = this.fMa;
            target.off = this.off;
        end
        
        function m = get.Ma(this)
            m = this.fMa;
        end
        
        function set.Ma(this, value)
            if ~isempty(find(isnan(value),1))
                error('Invalid coefficient values (NaN)');
            end
            this.machanged = true;
            this.fMa = value;
        end
        
        function set.off(this, value)
            if ~isempty(find(isnan(value),1))
                error('Invalid offset values (NaN)');
            end
            this.off = value;
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
            if ~isempty(this.off)
                copy.off = W' * this.off;
            end
            
        end
    end
    
end

