classdef BaseImplSolver < solvers.ode.BaseSolver
% BaseImplSolver: Base class for implicit ode solvers
%
% Some models like the models.pcd.PCDModel have a small CFL constant due to fast diffusion. Thus,
% implicit solvers are required to avoid long computation times.
%
% To not have to change the solve interface already present for explicit solvers, this class
% implements the solve method itself and uses a template method to perform the actual implicit
% solving.
%
% @author Daniel Wirtz @date 2011-04-19
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @new{0,3,dw,2011-04-19} Added this class to have a base for implicit solvers and wrap to the
% normal odefun interface from the explicit ones.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % A function handle to compute the core function's jacobian
        %
        % Optional, any implementations must also work with this property set to [].
        %
        % The function handles parameters must be:
        % - \c t The current time `t`
        % - \c x The current state space vector for time `t`
        % Returned must be the jacobian matrix `\frac{\partial f}{\partial x}(t,x)` at the point
        % `x`.
        %
        % @propclass{important} If available, supply a jacobian function evaluation handle to
        % improve speed and reliability of implicit solvers.
        %
        % See also: odeset
        JacFun = [];
    end
    
    methods
        
        function this = BaseImplSolver
            % Constructor. Sets the default name.
            this = this@solvers.ode.BaseSolver;
            
            this.Name = 'Base implicit solver';
            
            % Implicit solvers do not have a maximum time-step per default;
            % so overwrite the value here AFTER registration in BaseSolver
            % constructor.
            this.MaxStep = [];
            
            this.registerProps('JacFun');
        end
        
        function [t, x] = solve(this, odefun, t, x0)
            % Solves the ode specified by odefun implicitly.
            %
            % Parameters:
            % odefun: A function handle for the ODE's dynamic function `f(t,x)`
            % t: The times `t_i` on which to solve the ODE
            % x0: Initial condition vector `x_0(\mu)`
            %
            % Return values:
            % t: The desired computation times `t_0,\ldots,t_N` as row vector
            % x: The system's state `x_i` at time `t_i` as collection of column vectors
            implfun = @(t,x,xp)xp - odefun(t,x);
            
            opts = odeset;
            if ~isempty(this.MaxStep)
                opts = odeset(opts, 'MaxStep', this.MaxStep);
            end
            if ~isempty(this.InitialStep)
                opts = odeset(opts, 'InitialStep', this.InitialStep);
            end
            if ~isempty(this.JacFun)
                opts = odeset(opts, 'Jacobian', @this.FJAC);
            end
            
            [t,x] = this.implicit_solve(implfun, t, x0, odefun(0,x0), opts);
        end
        
        function set.JacFun(this, value)
            % Sets the jacobian function handle.
            %
            % Parameters:
            % value: a function handle that takes two arguments, `t` and `x`
            if ~isempty(value) && ~isa(value,'function_handle')
                error('JacFun must be a function handle');
            elseif ~isempty(value) && nargin(value) ~= 2
                error('JacFun must take exactly two arguments: t,x');
            end
            this.JacFun = value;
        end
    end
    
    methods(Access=private)
        function [dfdx,dfdxp] = FJAC(this, t, x, xp)%#ok
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
            
            % Implicit funcion is xp - odefun(t,x), so derivatives:
            dfdx = -this.JacFun(t, x);
            dfdxp = diag(ones(size(x, 1), 1));
        end
    end
    
    methods(Abstract)
        % Template method to perform actual implicit solving of the ODE.
        %
        % Parameters:
        % implfun: A handle to the implicit function `f` which describes
        % the ODE via `f(t,x,x') = 0` @type function_handle
        % t: The times `t_i` at which the ODE is to be computed @type rowvec
        % x0: Initial condition `x_0 = x(0)` @type colvec
        % xp0: Initial condition `x'_0 = x'(0)` @type colvec
        % opts: An odeset struct for additional options. @type struct
        %
        % Return values:
        % t: The times `t_i` @type rowvec
        % x: The solution `x(t_i)` of the ode at the times `t_i` @type matrix
        [t,x] = implicit_solve(this, implfun, t, x0, xp0, opts);
    end
    
end