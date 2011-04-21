classdef MLWrapper < solvers.ode.BaseSolver
    % Allows to wrap a MatLab ODE solver into the KerMor framework.
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
    
    methods
        
        function this = MLWrapper(solver)
            
            this = this@solvers.ode.BaseSolver;
            this.registerProps('MLSolver');
            
            if nargin > 0
                this.MLSolver = solver;
            end
        end
        
        function [t,y] = solve(this, odefun, t, x0)
            opts = [];
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
    
end

