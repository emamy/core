classdef qpIPOPT < solvers.BaseQPSolver
    % Quadratic program solver using IPOPT subroutines.
    %
    % @author Daniel Wirtz @date 27.10.2010
    
    properties
        % To provide user options accoring to the IPOPT MatlabInterface
        % specifications this property is a struct with fields that are
        % written to the "options" parameter struct of the ipopt main
        % routine. 
        % See @code help ipopt @endcode for options.
        UserOptions;
    end
    
    methods
        function this = qpIPOPT
            % Creates a qpIPOPT Wrapper for the Matlab-Interface.
            if isempty(which('ipopt'))
                error('ipopt installation not found. Added to path?');
            end
            this.Name = 'IPOPT QP Solver';
        end
        
        function copy = clone(this)
            copy = solvers.qpIPOPT;
            copy = clone@solvers.BaseQPSolver(this, copy);
            copy.UserOptions = this.UserOptions;
        end
    end
    
    methods(Access=protected)
        function [p,d,cflag,info] = internal_solve(this,Q,c,lb,ub,A,lbA,ubA,x0)
            
            m = size(Q,1);
            Q = sparse(Q);
            funcs = struct;
            %funcs.objective = @(beta).5*beta'*Q*beta + c'*beta;
            funcs.objective = @objective;
            
            if ~isempty(A)
                if ~issparse(A)
                    A = sparse(A);
                end
                funcs.gradient = @(beta)Q*beta + c + A';
                
                % A*beta = 0 constraints
                funcs.constraints = @(beta)A*beta;
                options.cl = lbA;
                options.cu = ubA;
                
                % Jacobian
                funcs.jacobian = @(beta)A;
                jstruct = sparse(ones(1,m));
                funcs.jacobianstructure = @()jstruct;
            else
                %funcs.gradient = @(beta)Q*beta + c;
                funcs.gradient = @gradient;
            end
            
            %funcs.hessian = @(beta,sigma,lambda)sparse(tril(sigma*Q));
            %hstruct = sparse(tril(ones(m)));
            funcs.hessian = @(beta,sigma,lambda)tril(sigma*Q);
            hstruct = sparse(double(tril(Q)>10*eps));
            funcs.hessianstructure = @()hstruct;
            
            % Bounds
            if ~isempty(lb)
                options.lb = lb;
            end
            if ~isempty(ub)
                options.ub = ub;
            end
            %options.bound_relax_factor = 0;
            
            %options.ipopt.print_level = KerMor.App.Verbose;
            %options.ipopt.output_file = 'out.txt';
            options.ipopt.tol         = 1e-7;
            
            if ~isempty(this.MaxIterations)
                options.ipopt.max_iter = this.MaxIterations;
            end
            
            %options.mu_strategy = 'monotone';
            %options.mu_oracle = 'probing';
            %options.corrector_type = 'affine';
            
            %options.mehrotra_algorithm = 'yes';
            %options.ipopt.linear_solver = 'ma27';
            %options.ipopt.linear_solver = 'ma57';
            %options.ipopt.linear_solver = 'mumps';
            
            %options.ipopt.hessian_approximation = 'limited-memory';
            %options.ipopt.derivative_test       = 'first-order';
            
            % Apply any user options and override internal settings if set.
            if ~isempty(this.UserOptions)
                if isstruct(this.UserOptions)
                    options = general.Utils.copyStructFields(this.UserOptions, options);
                else
                    warning('solvers:qpIPOPT','Property UserOptions must be a struct if set.');
                end
            end
            
            if isempty(x0)
                warning('solvers:qpIPOPT','No starting point given. Assuming Vector between bounds.');
                x0 = lb + (ub-lb)/2;
            end
            
            [p,info] = ipopt(x0,funcs,options);
            d = info.zl+info.zu;
            cflag = info.status == 0;
            info.Iterations = info.iter;
            
            function res = objective(beta)
                res = .5*beta'*Q*beta + c'*beta;
            end
            
            function res = gradient(beta)
                res = Q*beta + c;
            end
        end
    end    
end

