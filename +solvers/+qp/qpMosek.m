classdef qpMosek < solvers.qp.BaseQPSolver
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
        function this = qpMosek
            if isempty(which('mosekopt'))
                error('mosek installation not found. Added to path?');
            end
            this.Name = 'MOSEK 6.0 QP Solver';
        end
        
        function copy = clone(this)
            copy = solvers.qp.qpMosek;
            copy = clone@solvers.qp.BaseQPSolver(this, copy);
            copy.QuadProgOpts = this.QuadProgOpts;
        end
    end
    
    methods(Access=protected)
        function [p,d,cflag,info] = internal_solve(this,Q,c,lb,ub,A,lbA,ubA,x0)
            
            % Bounds
            Aeq = [];
            beq = [];
            Aineq = [];
            bineq = [];
            % Bounds
            if ~isempty(A)
                % Convert to equality and inequality bounds for this
                % interface (quadprog)
                eq = abs(lbA-ubA) < sqrt(eps);
                if any(eq)
                    Aeq = A(eq,:);
                    beq = lbA(eq);
                end
                if any(~eq)
                    Aineq = A(~eq,:);
                    bineq = ubA(~eq,:);
                end
            end
            
            % Set starting point if given
            s = [];
            if nargin == 9
                s = x0;
            end
            
            % Further options
            opts = this.QuadProgOpts;
            if ~isempty(this.MaxIterations)
                opts = optimset(opts,'MaxIter',this.MaxIterations);
            end
            
            % Solve QP H,f,A,b,B,c,l,u,x0,options
            [p, fval, exitflag, out, la] = mosek_quadprog(Q,c,Aineq,bineq,Aeq,beq,lb,ub,s,opts);
            
            % CARE! The -la.eqlin is subject to investigation and has to be
            % specified as soon as a suitable QP solver interface is
            % fixed/found
            cflag = exitflag > 0;
            d = [];
            if cflag > 0
                d = [la.lower + la.upper; -la.eqlin; la.ineqlin];
            end
            info = out;
            info.Iterations = out.iterations; %bummer!
            info.exitflag = exitflag;
        end
    end
    
end

