function ipopttest
Q = sparse([1 0 0; 0 1 0; 0 0 1]);
c = [1; 4; -1];
x0 = [1; 10; 0];
m = size(Q,1);

funcs = struct;
funcs.objective = @(beta).5*beta'*Q*beta + c'*beta;

%             if ~isempty(A)
%                 if ~issparse(A)
%                     A = sparse(A);
%                 end
%                 funcs.gradient = @(beta)Q*beta + c + A';
%
%                 % A*beta = 0 constraints
%                 funcs.constraints = @(beta)A*beta;
%                 options.cl = lbA;
%                 options.cu = ubA;
%
%                 % Jacobian
%                 funcs.jacobian = @(beta)A;
%                 jstruct = sparse(ones(1,m));
%                 funcs.jacobianstructure = @()jstruct;
%             else
funcs.gradient = @(beta)Q*beta + c;
%             end

funcs.hessian = @(beta,sigma,lambda)tril(sigma*Q);
hstruct = sparse(double(tril(Q)>10*eps));
funcs.hessianstructure = @()hstruct;

%options.bound_relax_factor = 0;

options.ipopt.print_level = 5;
%options.ipopt.output_file = 'out.txt';
options.ipopt.tol         = 1e-7;
options.ipopt.max_iter = 100;

%options.mu_strategy = 'monotone';
%options.mu_oracle = 'probing';
%options.corrector_type = 'affine';

%options.mehrotra_algorithm = 'yes';
%options.ipopt.linear_solver = 'ma27';
%options.ipopt.linear_solver = 'ma57';
%options.ipopt.linear_solver = 'mumps';

%options.ipopt.hessian_approximation = 'limited-memory';
%options.ipopt.derivative_test       = 'first-order';

[p,info] = ipopt(x0,funcs,options);
p
info

[p2,fval,exitflag,output,lambda] = quadprog(full(Q),c,[],[],[],[],[],[],x0,optimset('LargeScale','off'));
p2
exitflag
end

