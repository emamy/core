classdef MLode15i < solvers.MLWrapper & solvers.IImplSolver
% MLode15i: Wrapper for MatLab's ode15i builtin implicit solver
%
% @author Daniel Wirtz @date 2011-04-14
%
% @change{0,6,dw,2011-12-07} Changed the inheritance order, now inheriting
% from MLWrapper solver class and overriding the actual solver call
% (different arguments for ode15i)
%
% @change{0,5,dw,2011-10-16} Adopted to the new BaseSolver.RealTimeMode flag.
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @new{0,3,dw,2011-04-14} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetObservable)
        % Relative error tolerance for solver.
        %
        % Set to a lower value than the default 1e-3 of matlab ode solvers.
        %
        % See the RelTol description in the odeset documentation for the effects of this property.
        %
        % @propclass{critical} The correct relative tolerance guarantees the right number of
        % significant figures.
        %
        % @default 1e-4 @type double
        %
        % See also: AbsTol ode15i odeset
        RelTol = 1e-4;
        
        % Absolute error tolerance for solver.
        %
        % Set to a lower value than the default 1e-6 of matlab ode solvers.
        %
        % See the AbsTol description in the odeset documentation for the effects of this property.
        %
        % @propclass{critical} Absolute error tolerances require adoption to the current situation
        % and scale of computations.
        %
        % @default 1e-9 @type double
        %
        % See also: RelTol ode15i odeset
        AbsTol = 1e-9;
    end
    
    properties(Access=private)
        fED;
    end
        
    methods
        function this = MLode15i
            this = this@solvers.MLWrapper(@ode15i);
            this.Name = 'MatLab ode15i implicit solver wrapper';
            this.registerProps('RelTol','AbsTol');
            this.SolverType = solvers.SolverTypes.MLSolver;
            % "Disable" MaxStep DPCM warning as implicit solvers are stable
            this.MaxStep = [];
        end
    end
    
    methods(Access=protected)
        function varargout = solverCall(this, odefun, t, x0, opts)
            % Solves the ode specified by odefun implicitly.
            %
            % Parameters:
            % odefun: A function handle for the ODE's dynamic function `f(t,x)` @type
            % function_handle
            % t: The times `t_i` on which to solve the ODE @type rowvec<double>
            % x0: Initial condition vector `x_0(\mu)`
            % opts: Optional options struct for ODE settings obtained by \c odeset. @type
            % struct @default []
            %
            % Return values:
            % t: The desired computation times `t_0,\ldots,t_N` as row vector
            % x: The system's state `x_i` at time `t_i` as collection of column vectors
            opts = odeset(opts, 'RelTol', this.RelTol, 'AbsTol', this.AbsTol);
            
            %% Use properties from AJacobianSolver
            % Set Jacobian or Mass matrix
            if ~isempty(this.JacFun) || ~isempty(this.M)
                 opts = odeset(opts, 'Jacobian', @this.FJAC);
            end
            
            % Process any sparsity patterns
            JP = {[],[]};
            if ~isempty(this.JPattern)
                JP{1} = this.JPattern;
            end
            % Set df/dyp sparsity pattern, derived from mass matrix. Works
            % only (at least can be guaranteed) for constant mass matrices.
            if ~isempty(this.M)
                JP{2} = this.M.SparsityPattern;
            end
            opts = odeset(opts, 'JPattern', JP);
            
            %% Call implicit solver
            % implfun: A handle to the implicit function `f` which describes
            % the ODE via `f(t,x,x') = 0`
            xp0 = odefun(0,x0);
            if ~isempty(this.M)
                implfun = @(t,x,xp)this.M.evaluate(t)*xp - odefun(t,x);
                xp0 = this.M.evaluate(0)\xp0;
            else
                implfun = @(t,x,xp)xp - odefun(t,x);
            end
            % Call ode15i solver
            [varargout{1:nargout}] = this.MLSolver(implfun, t, x0, xp0, opts);
        end
    end
    
    methods(Access=private)
        function [dfdx, dfdxp] = FJAC(this, t, x, xp)%#ok
            % Internal implementation of the FJAC function utilized by i.e. MatLab's builtin
            % solvers.
            %
            % Parameters:
            % t: The time `t`
            % x: The current solution vector `x(t)`
            % xp: The current derivative solution vector `x'(t)`
            %
            % Return values:
            % dfdx: The jacobian of the implicit function `f(t,x,x')` with respect to `x`
            % dfdxp: The jacobian of the implicit function `f(t,x,x')` with respect to `x'`
            % This value corresponds to mass matrices `M` for systems of the
            % type `Mx'(t) = f(x(t),t) \ldots`
            
            % Implicit funcion is M*xp - odefun(t,x), so derivatives:
            dfdx = [];
            if ~isempty(this.JacFun)
                dfdx = -this.JacFun(t, x);
            end
            dfdxp = [];
            if ~isempty(this.M)
                dfdxp = this.M.evaluate(t);
            end
        end
    end
    
end