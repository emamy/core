classdef BaseQPSolver < KerMorObject & ICloneable
    % Base class for any KerMor QP-Solver or QP-Wrapper
    %
    % Solvers solve quadratic programming problems of the form
    % `\frac{1}{2}\beta Q\beta + c^t\beta` for `Q` positive definite matrices `Q`.
    %
    % @author Daniel Wirtz @date 27.10.2010
    %
    % @change{0,4,dw,2011-05-05} Included this class into @ref propclasses. Also set a default value
    % for solvers.qp.BaseQPSolver.MaxIterations and a suggestion (For Verbose > 0) if it might be too small.
    %
    % @todo write tests for solvers incorporating different constraint
    % settings (equality, inequality, upper, lower etc)
    
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
        
        % The name of the Solver
        %
        % @propclass{optional} The qp solver's name.
        %
        % @type char @default BaseQPSolver
        Name = 'BaseQPSolver';
    end
    
    methods
        
        function this = BaseQPSolver
            this = this@KerMorObject;
            
            this.registerProps('MaxIterations','Name');
        end
        
        function target = clone(this, target)
            if nargin == 1
                error('Cannot call clone without subclass argument.');
            end
            target.MaxIterations = this.MaxIterations;
            target.Name = this.Name;
        end
        
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
                    fprintf('%s: MaxIterations of %d smaller than recommended value %d.\n',this.Name,this.MaxIterations,it);
                end
            end
            
            if nargin < 9
                x0 = [];
            end
            
            s = tic;
            % Call internal solver algorithm
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
                m = MException('KerMor:solvers:qp:notconverged','Solver ''%s'' did not converge.',this.Name);
                m.throw;
            end
            
            % Verbose output
            if KerMor.App.Verbose > 3
                fprintf('QP-Solver ''%s'' finished after %d/%d Iterations and %f seconds. Info:\n',this.Name,info.Iterations,this.MaxIterations, info.CompTime);
                disp(info);
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Template method for specific quadratic program solution
        % algorithms.
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
        % cflag: Boolean. Set to true if satisfying result was obtained
        % info: Struct that contains algorithm-specific output information.
        % Must contain a field "Iterations"
        [p,d,cflag,info] = internal_solve(this,Q,c,lb,ub,A,Alb,Aub,x0);
    end
    
    methods(Static)
        function res = test_QPSolversNoBounds
            
%             testsize = 50;
%             k = kernels.GaussKernel(10);
%             x = linspace(-5,5,testsize);
            
            Q = [1 0 0; 0 1 0; 0 0 1];
            c = [1; 4; -1];
            x0 = [1; 50; 0];
            
            qp{1} = solvers.qp.qpMatlab;
            qp{2} = solvers.qp.qpOASES;
            qp{3} = solvers.qp.qpMosek;
%             ip = solvers.qp.qpIPOPT;
%             ip.UserOptions.ipopt.linear_solver = 'ma27';
%             qp{4} = ip;
%             ip = solvers.qp.qpIPOPT;
%             ip.UserOptions.ipopt.linear_solver = 'ma57';
%             qp{5} = ip;
%             ip = solvers.qp.qpIPOPT;
%             ip.UserOptions.ipopt.linear_solver = 'mumps';
%             qp{6} = ip;
            
            for idx = 1:length(qp)
                %qp{idx}.MaxIterations = 50;
                [p(:,idx),d(:,idx),info{idx}] = qp{idx}.solve(Q,c,[],[],[],[],[],x0);
                info{idx}
            end
            % Success only if all results are very similar
            res = norm(repmat(p(:,1),1,size(p,2))-p) < sqrt(eps);
        end
    end
    
end

