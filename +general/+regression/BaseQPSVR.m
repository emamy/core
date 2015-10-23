classdef BaseQPSVR < general.regression.BaseScalarSVR
% BaseQPSVR: SVR variant that is solved using quadratic programs.
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-01-23
%
% @new{0,7,dw,2013-01-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetObservable)
        % The maximum number of iterations.
        %
        % If set to [], the number is automatically set to
        % <code>100*(size(Q,1)+size(A,1));</code>
        %
        % @propclass{alglimit} Works as an execution limit when the qp does not converge.
        %
        % @default 5000 @type integer
        MaxIterations = 5000;
        
        % Options for quadprog-solver
        %
        % @propclass{important} The flags for the matlab builtin quadprog solver.
        %
        % See also: quadprog
        QuadProgOpts = optimset('Display','off','Algorithm','interior-point-convex');
        %'LargeScale','off', interior-point-convex, trust-region-reflective
    end
    
    methods
        function this = BaseQPSVR
            this = this@general.regression.BaseScalarSVR;
            
            this.registerProps('MaxIterations','QuadProgOpts');
        end
        
        function copy = clone(this, copy)
            if nargin == 1
                error('Cannot call clone without subclass argument.');
            end
            copy = clone@general.regression.BaseScalarSVR(this, copy);
            copy.MaxIterations = this.MaxIterations;
            copy.QuadProgOpts = this.QuadProgOpts;
        end
    end
    
    methods(Access=protected)
        function [p,d,info] = solve(this,Q,c,lb,ub,A,Alb,Aub,x0)
            % Solves the given quadratic problem `\frac{1}{2}\beta Q\beta +
            % c^t\beta` according to the subclasses' algorithm.
            %
            % Throws an error if the solver does not come to a solution
            % satisfying the requirements (MaxIterations, Tolerance etc
            % depending on subclasses)
            %
            % Parameters:
            % Q: The quadratic pos. def. matrix Q @type matrix
            % c: The linear part `c` @type colvec
            % lb: The lower bounds of `\beta` @type colvec
            % ub: The upper bounds of `\beta` @type colvec
            % A: The equality constraint matrix `A\beta=0` @type matrix
            % Alb: The lower bound matrix `lb < A_{lb}\beta` @type matrix
            % Aub: The upper bound matrix `A_{ub}\beta < ub` @type matrix
            % x0: The initial value for `\beta`
            %
            % Return values:
            % p: The primary variable `\beta`
            % d: The dual variables, i.e. the LaGrange multipliers for the
            % constraints. Correspond to a column vector composed as
            % [lambda_bounds; lambda_equality; lambda_inequality]
            % info: Contains algorithm-specific return values, at minimum
            % the self-explaining fields 
            % # Iterations
            % # CompTime
            % # NumVariables
            % # NumConstraints
            %
            % @throws KerMor:solvers:qp:notconverged Thrown if the QP solver did not converge.
            
            % Check for MaxIterations
            if KerMor.App.Verbose > 0
                it = 100*(size(Q,1)+size(A,1));
                if KerMor.App.Verbose > 3 && this.MaxIterations < it
                    fprintf('BaseQPSVR: MaxIterations of %d smaller than recommended value %d.\n',this.MaxIterations,it);
                end
            end
            
            if nargin < 9
                x0 = [];
            end
            
            s = tic;
            % Call internal solver algorithm (currently: only matlab's quadprog)
            [p,d,cflag,info] = this.internal_solve(Q,c,lb,ub,A,Alb,Aub,x0);
            
            % Adds a CompTime field
            info.CompTime = toc(s);
            info.NumVariables = size(Q,1);
            info.NumConstraints = size(A,1);
            
            % Evaluate
            if ~cflag
                disp(info);
                if isfield(info,'message')
                    fprintf('Message: %s\n',info.message);
                end
                m = MException('BaseQPSVR:notconverged','quadprog solver did not converge.');
                m.throw;
            end
            
            % Verbose output
            if KerMor.App.Verbose > 2
                fprintf('QP-SVR finished after %d/%d Iterations and %f seconds. Info:\n',info.Iterations,this.MaxIterations, info.CompTime);
                disp(info);
            end
        end
    end
    
    methods(Access=private)
        function [p,d,cflag,info] = internal_solve(this,Q,c,lb,ub,A,lbA,ubA,x0)
            % Solves the QP using MatLabs ''quadprog'' routine
            
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
            [p, ~, exitflag, out, la] = quadprog(prog);
            
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