classdef qpMatlab < solvers.qp.BaseQPSolver
    % Quadratic program solver using matlab's quadprog subroutine.
    %
    % The specific behaviour of quadprog can be steered by using the 
    % QuadProgOpts property.
    %
    % @author Daniel Wirtz @date 27.10.2010
    
    properties
        % Options for quadprog-solver
        QuadProgOpts = optimset('LargeScale','off','Display','off');
    end
    
    methods
        function this = qpMatlab
            this.Name = 'Matlab native quadprog solver';
        end
        
        function copy = clone(this)
            copy = solvers.qp.qpMatlab;
            copy = clone@solvers.qp.BaseQPSolver(this, copy);
            copy.QuadProgOpts = this.QuadProgOpts;
        end
    end
    methods(Access=protected)
        function [p,d,cflag,info] = internal_solve(this,Q,c,lb,ub,A,lbA,ubA,x0)
            
            % Problem setup
            prog.H = Q;
            prog.f = c;
            
            % Bounds
            if ~isempty(A)
                % Convert to equality and inequality bounds for this
                % interface (quadprog)
                eq = abs(lbA-ubA) < sqrt(eps);
                if any(eq)
                    prog.Aeq = A(eq,:);
                    prog.beq = lbA(eq);
                end
                if any(~eq)
                    prog.Aineq = A(~eq,:);
                    prog.bineq = ubA(~eq,:);
                end
            end
            prog.lb = lb;
            prog.ub = ub;
            
             % Set starting point if given
            if nargin == 9
                prog.x0 = x0;
            end
            
            % Further options
            opts = this.QuadProgOpts;
            if ~isempty(this.MaxIterations)
               opts = optimset(opts,'MaxIter',this.MaxIterations);
            end
            prog.options = opts;
            prog.solver = 'quadprog';
            
            % Solve QP
            [p, fval, exitflag, out, la] = quadprog(prog);
            
            % CARE! The -la.eqlin is subject to investigation and has to be
            % specified as soon as a suitable QP solver interface is
            % fixed/found
            d = [la.lower + la.upper; -la.eqlin; la.ineqlin];
            cflag = exitflag == 1;
            info = out;
            info.exitflag = exitflag;
            info.Iterations = out.iterations; %bummer!
        end
    end
    
end

