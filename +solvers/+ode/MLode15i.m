classdef MLode15i < solvers.ode.BaseImplSolver
% MLode15i: Wrapper for MatLab's ode15i builtin implicit solver
%
% @author Daniel Wirtz @date 2011-04-14
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
            this = this@solvers.ode.BaseImplSolver;
            this.Name = 'MatLab ode15i implicit solver wrapper';
            this.registerProps('RelTol','AbsTol');
        end
        
        function [t,x] = implicit_solve(this, implfun, t, x0, x0p, opts)
            % Solves the ode implicitly.
            %
            % Parameters:
            % implfun: A handle to the implicit function `f` which describes
            % the ODE via `f(t,x,x') = 0` @type function_handle
            % t: The times `t_i` at which the ODE is to be computed @type rowvec
            % x0: Initial condition `x_0 = x(0)` @type colvec
            % x0p: Initial condition `x'_0 = x'(0)` @type colvec
            % opts: An odeset struct for additional options. @type struct
            if nargin < 6
                opts = odeset;
            end
            opts = odeset(opts, 'RelTol', this.RelTol, 'AbsTol', this.AbsTol);
            opts = odeset(opts, 'OutputFcn',@this.ODEOutputFcn);
            this.fED = solvers.ode.SolverEventData;
            [t,x] = ode15i(implfun, t, x0, x0p, opts);
            x = x';
            t = t';
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