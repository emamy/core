classdef MLWrapper < solvers.BaseSolver
    % Allows to wrap a MatLab ODE solver into the KerMor framework.
    
    properties
        % The wrapped Matlab-Solver
        % Function handle.
        %
        % Default: ode23
        % See also: ode23 ode45
        MLSolver = @ode23;
    end
    
    methods
        
        function this = MLWrapper(solver)
            % Creates the wrapper class taking 
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

