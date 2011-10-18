classdef MLWrapper < solvers.ode.BaseSolver
    % Allows to wrap a MatLab ODE solver into the KerMor framework.
    %
    % @author Daniel Wirtz @date 2010-08-09
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
    end
    
    properties(Access=private)
        % Event data for StepPerformed
        fED;
    end
    
    methods
        
        function this = MLWrapper(solver)
            this = this@solvers.ode.BaseSolver;
            this.registerProps('MLSolver');
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
            this.fED = solvers.ode.SolverEventData;
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

