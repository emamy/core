classdef BaseQPSolver < ICloneable
    % Base class for any KerMor QP-Solver or QP-Wrapper
    %
    % Solvers solve quadratic programming problems of the form
    % `\frac{1}{2}\beta Q\beta + c^t\beta` for `Q` positive definite
    % @docupdate
    %
    % @todo write tests for solvers incorporating different constraint
    % settings (equality, inequality, upper, lower etc)
    %
    %
    % @author Daniel Wirtz @date 27.10.2010
    
    properties
        % The maximum number of iterations.
        % If set to [], the number is automatically set to
        % <code>100*(size(Q,1)+size(A,1));</code>
        MaxIterations;
        
        % The name of the Solver
        Name = 'Not specified - BaseQPSolver';
    end
    
    methods
        
        function target = clone(this, target)
            if nargin == 1
                error('Cannot call clone without subclass argument.');
            end
            target.MaxIterations = this.MaxIterations;
            target.Name = this.Name;
        end
        
        function [p,d,info] = solve(this,Q,c,lb,ub,A,Alb,Aub,x0)
            % Solves the given quadratic problem according to the
            % subclasses' algorithm.
            %
            % Throws an error if the solver does not come to a solution
            % satisfying the requirements (MaxIterations, Tolerance etc
            % depending on subclasses)
            %
            % Parameters:
            % Q:
            % c:
            %
            % Return values:
            % p:
            % d: The dual variables, i.e. the LaGrange multipliers for the
            % constraints. Correspond to a column vector composed as
            % [lambda_bounds; lambda_equality; lambda_inequality]
            % info: Contains algorithm-specific return values, at minimum
            % the self-explaining fields 
            % # Iterations
            % # CompTime
            % # NumVariables
            % # NumConstraints
            
            % Check for MaxIterations
            it = this.MaxIterations;
            if isempty(it)
                it = 100*(size(Q,1)+size(A,1));
                fprintf('%s: No MaxIterations set, assuming MaxIterations=%d\n',this.Name,it);
                this.MaxIterations = it;
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
                error('Solver ''%s'' did not converge.',this.Name);
            end
            
            % Verbose output
            if KerMor.Instance.Verbose > 3
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
        % x0: If none is given, [] is passed.
        %
        % Return values:
        % cflag: Boolean. Set to true if satisfying result was obtained
        % info: Struct that contains algorithm-specific output information.
        % Must contain a field "Iterations"
        [p,d,cflag,info] = internal_solve(Q,c,lb,ub,A,Alb,Aub,x0);
    end
    
    methods(Static)
        function res = test_QPSolversNoBounds
            
%             testsize = 50;
%             k = kernels.GaussKernel(10);
%             x = linspace(-5,5,testsize);
            
            Q = [1 0 0; 0 1 0; 0 0 1];
            c = [1; 4; -1];
            x0 = [1; 50; 0];
            
            qp{1} = solvers.qpMatlab;
            qp{2} = solvers.qpOASES;
            qp{3} = solvers.qpMosek;
            ip = solvers.qpIPOPT;
            ip.UserOptions.ipopt.linear_solver = 'ma27';
            qp{4} = ip;
            ip = solvers.qpIPOPT;
            ip.UserOptions.ipopt.linear_solver = 'ma57';
            qp{5} = ip;
            ip = solvers.qpIPOPT;
            ip.UserOptions.ipopt.linear_solver = 'mumps';
            qp{6} = ip;
            
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

