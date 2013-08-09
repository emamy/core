classdef MLWrapper < solvers.BaseSolver & solvers.AJacobianSolver
% Allows to wrap a MatLab ODE solver into the KerMor framework.
%
% @author Daniel Wirtz @date 2010-08-09
%
% @new{0,7,dw,2013-05-16} Now inheriting from the new
% solvers.AJacobianSolver abstract class to all MatLab wrapper solvers
% (e.g. ode15s also benefits from the Jacobian)
%
% @new{0,6,dw,2011-12-07} Added a new field odeopts which corresponds to
% the odeset options struct of MatLab. Any values set will be passed to the
% internal ode23 etc solver.
%
% @change{0,5,dw,2011-10-16} Adopted to the new BaseSolver.RealTimeMode flag.
%
% @change{0,5,dw,2011-10-15} Moved the creation of the SolverEventData
% into the solve function as creation in the constructor seems to crash
% Matlab versions prior to the 2011a which was used to program this
% functionality in the first place.
%
% @change{0,5,dw,2011-09-29} Added callback for StepPerformed to enable
% "real time" plotting.
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The wrapped Matlab-Solver
        %
        % Has to be a function handle to one of Matlab's implemented ODE
        % solvers.
        %
        % @propclass{critical} The correct underlying MatLab builtin solver can make the difference.
        %
        % @default ode23 @type function_handle
        %
        % See also: ode23 ode45
        MLSolver = @ode23;
        
        % Additional ODE options.
        %
        % Change this value via the odeset function. This struct will be
        % passed as a whole to the wrapped MatLab solver.
        %
        % @type struct @default odeset
        %
        % See also: odeset
        odeopts;
    end
    
    properties(Access=private)
        % Event data for StepPerformed
        fED;
    end
    
    methods
        
        function this = MLWrapper(solver)
            this = this@solvers.BaseSolver;
            this.registerProps('MLSolver');
            this.SolverType = solvers.SolverTypes.MLSolver;
            if nargin > 0
                this.MLSolver = solver;
            end
            this.odeopts = odeset;
        end
        
        function [t,y] = solve(this, odefun, t, x0)
            opts = this.odeopts;
            if ~isempty(this.MaxStep)
                opts = odeset(opts, 'MaxStep', this.MaxStep);
            end
            if ~isempty(this.InitialStep)
                opts = odeset(opts, 'InitialStep', this.InitialStep);
            end
            
            % Use AJacobianSolver properties
            opts = odeset(opts, 'Jacobian', this.JacFun);
            JP = [];
            % Matlab solvers only use patterns if the jacobian function is
            % not directly given
            if isempty(this.JacFun)
                JP = this.JPattern;
            end
            opts = odeset(opts, 'JPattern', JP);
                
            % Pass Mass Matrix to solver (only applicable for non-ode15i
            % solvers, the latter one makes use of M in a different way,
            % see the class solvers.MLode15i)
            M = []; yp0 = []; MS = [];
            if ~isempty(this.M)
                if ~this.M.TimeDependent
                    M = this.M.evaluate(0);
                    opts = odeset(opts,'MassConstant','true');
                else
                    M = @(t)this.M.evaluate(t);
                end
                % Compute initial slope
                yp0 = this.M.evaluate(0)\odefun(0,x0);
                MS = 'none';
            end
            opts = odeset(opts,'Mass',M,'MStateDependence',MS,...
                    'InitialSlope',yp0);
            
            if this.RealTimeMode
                opts = odeset(opts,'OutputFcn',@this.ODEOutputFcn);
                % Bug in Matlab 2009a: direct assignment crashes Matlab!
                % Seems also not to work if created within the constructor.
                ed = solvers.SolverEventData;
                this.fED = ed;
                this.solverCall(odefun, t, x0, opts);
                t = []; y = [];
            else
                [t,y] = this.solverCall(odefun, t, x0, opts);
                y = y';
                t = t';
            end
        end
        
        function set.MLSolver(this, value)
            if isa(value,'function_handle')
                this.MLSolver = value;
                this.Name = sprintf('Matlab Solver wrapper using %s',func2str(value));
            else
                error('Invalid function handle!');
            end
        end
    end
    
    methods(Access=protected)
        function varargout = solverCall(this, odefun, t, x0, opts)
            % Default solver call for all builtin ode solvers except
            % ode15i. This method gets overridden in MLode15i.
            %
            % Parameters:
            % odefun: A function handle for the ODE's dynamic function `f(t,x)` @type
            % function_handle
            % t: The times `t_i` on which to solve the ODE @type rowvec<double>
            % x0: Initial condition vector `x_0(\mu)`
            % opts: Optional options struct for ODE settings obtained by \c odeset. @type
            % struct @default []
            [varargout{1:nargout}] = this.MLSolver(odefun, t, x0, opts);
        end
    end
    
    methods(Access=protected, Sealed)
        function status = ODEOutputFcn(this, t, y, flag)
            % Wraps the OutputFcn of the Matlab ODE solver into
            % the StepPerformed event
            %
            % Parameters:
            % t: The current time `t`
            % y: The system's output `y(t)`
            % flag: The flag passed from the ODE solver as argument of the
            % 'OutputFcn' setting in odeset.
            %
            % See also: odeset
            if ~strcmp(flag,'init')
                % For some reason the t and y args have more than one
                % entry, so loop over all of them.
                for idx=1:length(t)
                    this.fED.Times = t(idx); % when flag==init the t var is larger than one
                    this.fED.States = y(:,idx);
                    this.notify('StepPerformed',this.fED);
                end
            end
            status = 0;
        end
    end
    
end

