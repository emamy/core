classdef ParamTimeKernelCoreFun < kernels.ParamTimeKernelExpansion & dscomponents.ACoreFun
% ParamTimeKernelCoreFun: Dynamical system core function based on an parametric time dependent
% kernel expansion.
%
% Convenience class that wraps a given kernel expansion into the ACoreFun interface.
% The main methods 'clone' and 'evaluate' are redefined with explicit choice of which one to
% inherit.
%
%
%
% @author Daniel Wirtz @date 2011-07-07
%
% @new{0,5,dw,2011-07-07} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = ParamTimeKernelCoreFun
            this = this@kernels.ParamTimeKernelExpansion;
            this = this@dscomponents.ACoreFun;
            
            this.CustomProjection = true;
            this.MultiArgumentEvaluations = true;
            this.TimeDependent = false;
            this.CustomJacobian = true;
            this.addlistener('TimeKernel','PostSet',@this.TimeKernelPostSet);
        end
        
        function projected = project(this, V, W)
            % Call superclass method
            projected = this.clone; 
            
            projected = project@dscomponents.ACoreFun(this, V, W, projected);
            % For rotation invariant kernel expansions the snapshots can be
            % transferred into the subspace without loss.
            if this.RotationInvariant
                projected.Centers.xi = W' * this.Centers.xi;
            end
            projected.Ma = W'*this.Ma;
        end
        
        function fx = evaluate(this, x, t, mu)
            % From both inherited evaluate functions take the kernel evaluation, augmented by the
            % `V` multiplication if projection is used but the kernel is not rotation invariant.
%             V = 1;
%             if ~this.RotationInvariant && ~isempty(this.V)
%                 V = this.V;
%             end
            fx = evaluate@kernels.ParamTimeKernelExpansion(this, x, t, mu);
        end
        
        function phi = getKernelVector(this, x, t, mu)
            V = 1;
            if ~this.RotationInvariant && ~isempty(this.V)
                V = this.V;
            end
            phi = getKernelVector@kernels.ParamTimeKernelExpansion(this, V*x, t, mu);
        end
        
        function copy = clone(this, copy)
            if nargin < 2
                copy = dscomponents.ParamTimeKernelCoreFun;
            end
            copy = clone@kernels.ParamTimeKernelExpansion(this, copy);
            copy = clone@dscomponents.ACoreFun(this, copy);
        end
        
        function J = getStateJacobian(this, x, t, mu)
            % Implement explicitly as both ACoreFun and KernelExpansion
            % provide getStateJacobian methods.
            J = getStateJacobian@kernels.ParamTimeKernelExpansion(this, x, t, mu);
        end
        
        function y = evaluateCoreFun(this)%#ok
            % Noting to do here, evaluate is implemented directly. This method will never be called.
        end
    end
    
    methods(Access=private)
        function TimeKernelPostSet(this, ~, ~)
            this.TimeDependent = ~isa(this.TimeKernel,'kernels.NoKernel');
        end
    end
    
    methods(Static,Access=protected)
        function this = loadobj(this)
            this = loadobj@DPCMObject(this);
            % Register listener for TimeDependent changed
            this.addlistener('TimeKernel','PostSet',@this.TimeKernelPostSet);
        end
    end
    
end