classdef ParamTimeKernelCoreFun < dscomponents.ACoreFun & dscomponents.IGlobalLipschitz
% ParamTimeKernelCoreFun: Dynamical system core function which evaluates a
% contained kernel expansion, either parametric or plain state-dependence.
%
% Convenience class that wraps a given kernel expansion into the ACoreFun
% interface.
% The main methods 'clone' and 'evaluate' are redefined with explicit
% choice of which one to inherit.
%
% @author Daniel Wirtz @date 2011-07-07
%
% @change{0,7,dw,2014-01-15} Removed the inheritance from
% kernels.ParamTimeKernelExpansion and included an instance instead. (Favor
% inclusion over inheritance, Warwick said!)
%
% @new{0,5,dw,2011-07-07} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    properties(SetObservable)
        % The inner kernel expansion which is evaluated as core function.
        %
        % @type kernels.KernelExpansion @default
        % kernels.ParamTimeKernelExpansion
        Expansion;
    end
    
    properties(Access=private, Transient)
        % Listener for centers change
        cl;
        % Listener for coefficients change
        mal;
    end
    
    methods
        function this = ParamTimeKernelCoreFun(sys)
            this = this@dscomponents.ACoreFun(sys);
            
            this.CustomProjection = true;
            this.TimeDependent = false;
            this.addlistener('Expansion','PostSet',@this.ExpansionPostSet);
            % addDimListeners is triggered in Expansion post set!
            this.Expansion = kernels.ParamTimeKernelExpansion;
        end
        
        function projected = project(this, V, W)
            % Call superclass method
            projected = this.clone;
            projected = project@dscomponents.ACoreFun(this, V, W, projected);
            % For rotation invariant kernel expansions the snapshots can be
            % transferred into the subspace without loss.
            lkexp = this.Expansion;
            kexp = lkexp.clone;
            k = kexp.Kernel;
            if k.IsRBF
                kexp.Centers.xi = W' * lkexp.Centers.xi;
                PG = V'*(k.G*V);
                if all(all(abs(PG - eye(size(PG,1))) < 100*eps))
                    PG = 1;
                end
                kexp.Kernel.G = PG;
            elseif k.IsScProd
                kexp.Centers.xi = V' * (k.G * lkexp.Centers.xi);
                kexp.Kernel.G = 1;
            end
            kexp.Ma = W'*lkexp.Ma;
            projected.Expansion = kexp;
        end
        
        function fx = evaluate(this, x, t)
            % Evaluates this CoreFun
            %
            % Parameters:
            % x: The current state space position @type colvec<double>
            % varargin: For ParamTimeKernelExpansions, additionally `t` and
            % `\mu` can be provided.
            % t: The time `t` @type double
            % mu: The parameter `\mu` @type colvec<double>
            %
            % Return values:
            % fx: The evaluation of the kernel expansion @type
            % matrix<double>
            if nargin < 3
                t = 0;
            end
            fx = this.Expansion.evaluate(x, t, this.mu);
        end
        
        function fx = evaluateMulti(this, x, varargin)
            % Evaluates this CoreFun for multiple values
            %
            % Parameters:
            % x: The state space positions @type matrix<double>
            % varargin: For ParamTimeKernelExpansions, additionally `t` and
            % `\mu` can be provided.
            % t: The times `t` @type rowvec<double>
            % mu: The parameters `\mu` @type matrix<double>
            %
            % Return values:
            % fx: The evaluation of the kernel expansion on all data @type
            % matrix<double>
            fx = this.Expansion.evaluate(x, varargin{:});
        end
        
        function L = getGlobalLipschitz(this, t, mu)
            L = this.Expansion.getGlobalLipschitz(t, mu);
        end
               
        function copy = clone(this, copy)
            if nargin < 2
                copy = dscomponents.ParamTimeKernelCoreFun(this.System);
            end
            copy.Expansion = this.Expansion.clone;
            copy = clone@dscomponents.ACoreFun(this, copy);
        end
        
        function J = getStateJacobian(this, x, varargin)
            % Implement explicitly as both ACoreFun and KernelExpansion
            % provide getStateJacobian methods.
            %
            % Parameters:
            % x: The current state space position @type colvec<double>
            % varargin: For ParamTimeKernelExpansions, additionally `t` and
            % `\mu` can be provided.
            % t: The time `t` @type double
            % mu: The parameter `\mu` @type colvec<double>
            %
            % Return values:
            % J: The state jacobian @type matrix<double>
            J = this.Expansion.getStateJacobian(x, varargin{:});
        end
        
        function y = evaluateCoreFun(this)%#ok
            error('This should never be called. "evaluate" is implemented directly.');
            % Noting to do here, evaluate is implemented directly. This method will never be called.
        end
    end
    
    methods(Access=private)
        function ExpansionPostSet(this, ~, ~)
            this.fDim = 0;
            this.xDim = 0;
            if ~isempty(this.Expansion)
                this.TimeDependent = ...
                    isa(this.Expansion, 'kernels.ParamTimeKernelExpansion') ...
                    && ~isa(this.Expansion.TimeKernel,'kernels.NoKernel');
                this.fDim = size(this.Expansion.Ma,1);
                this.xDim = size(this.Expansion.Centers.xi,1);
            end
            this.addDimListeners;
        end
        
        function addDimListeners(this)
            delete(this.cl);
            delete(this.mal);
            if ~isempty(this.Expansion)
                this.cl = this.Expansion.addlistener('Centers','PostSet',@this.CentersPostSet);
                this.mal = this.Expansion.addlistener('Ma','PostSet',@this.MaPostSet);
            end
        end
        
        function CentersPostSet(this, ~, ~)
            this.xDim = size(this.Expansion.Centers.xi,1);
        end
        
        function MaPostSet(this, ~, ~)
            this.fDim = size(this.Expansion.Ma,1);
        end
    end
    
    methods(Static,Access=protected)
        function this = loadobj(this)
            this = loadobj@DPCMObject(this);
            % Register listener for TimeDependent changed
            this.addlistener('Expansion','PostSet',@this.ExpansionPostSet);
            this.addDimListeners;
        end
    end
    
end