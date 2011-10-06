classdef MLWrapper < solvers.ode.BaseSolver
    % Allows to wrap a MatLab ODE solver into the KerMor framework.
    %
    % @change{0,5,dw,2011-09-29} Added callback for StepPerformed to enable
    % "real time" plotting.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    
    properties(SetObservable)
        % The wrapped Matlab-Solver
        % Function handle.
        %
        % @propclass{critical} The correct underlying MatLab builtin solver can make the difference.
        %
        % @default ode23
        %
        % See also: ode23 ode45
        MLSolver = @ode23;
    end
    
    properties(Access=private)
        % Event data for StepPerformed
        fED;
    end
    
    methods
        
        function this = MLWrapper(solver)
            
            this = this@solvers.ode.BaseSolver;
            this.registerProps('MLSolver');
            this.fED = solvers.ode.SolverEventData;
            if nargin > 0
                this.MLSolver = solver;
            end
        end
        
        function [t,y] = solve(this, odefun, t, x0)
            opts = odeset('OutputFcn',@this.ODEOutputFcn);
            if ~isempty(this.MaxStep)
                opts = odeset(opts, 'MaxStep', this.MaxStep);
            end
            if ~isempty(this.InitialStep)
                opts = odeset(opts, 'InitialStep', this.InitialStep);
            end
            [t,y] = this.MLSolver(odefun, t, x0, opts);
            y = y';
            t = t';
        end
        
        function set.MLSolver(this, value)
            if isa(value,'function_handle')
                this.MLSolver = value;
                this.Name = sprintf('Matlab Solver wrapper using %s',func2str(value));%#ok
            else
                error('Invalid function handle!');
            end
        end
    end
    
    methods(Access=protected, Sealed)
        function status = ODEOutputFcn(this, t, y, flag)
            % Wraps the OutputFcn of the Matlab ODE solver into
            % the StepPerformed event
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

