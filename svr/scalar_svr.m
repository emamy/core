function [ai,b,svidx] = scalar_svr(fxi,kernelmatrix,eps,C)
%SCALAR_SVR Performs scalar support vector regression
%   TODO

if nargin < 4
    C = 1;
end

% Maximum number of iterations for automatic b-computation in case all
% alpha_i's are bounded by C/m
maxIt = 100;

% Compile quadratic program
% Total number of samples
m = size(kernelmatrix,1);
% Ensure fxi is a column vector
fxi = reshape(fxi,m,[]);

% Storage: alpha(1..m) = alpha_i, alpha(m+1..2m) = alpha_i*
% T performs alpha_i* - alpha_i each
T = [-diag(ones(1,m)) diag(ones(1,m))];

%% Problem setup
prog.H = T'*kernelmatrix*T;
prog.f = eps*ones(2*m,1) - T'*fxi;

prog.Aeq = [ones(1,m) -ones(1,m)];
prog.beq = 0;

prog.lb = zeros(2*m,1);
prog.ub = ones(2*m,1)*(C/m);

prog.Aineq = [];
prog.bineq = [];
%prog.x0 = rand(2*m,1)*(C/m);

prog.solver = 'quadprog';
prog.options = optimset('LargeScale','off','MaxIter',300,'Display','off');

%% Solve quadratic problem
ai = quadprog(prog);

%% Extract support vectors
% reduce ai from 2m to m vector
% follow alpha_i* - alpha_i
ai = T*ai;

% Find support vectors
svidx = find(abs(ai) >= 1e-10);

if isempty(svidx)
    error('No support vectors found. Problem unsolvable with current config?');
end
% Reduce set to support vectors


% check if b can be computed correctly
if all(abs(ai(svidx)) - C/m < 1e-10)
    b = 0;
    
    % Save to-leave-out vectors
    skipped = setdiff(1:m,svidx);
    fxiskip = fxi(skipped)';
    
    ai = ai(svidx);
    
    oldb = Inf;
    cnt = 0;
    while abs(oldb-b) > 1e-5 && cnt < maxIt
        % Compute approx fx values at xskipped
        fsvrskip = ai'*kernelmatrix(svidx,skipped) + b;
        % Comput max and min differences
        up  = max(fxiskip - fsvrskip + eps);
        low = min(fxiskip - fsvrskip - eps);
        
        oldb = b;
        b = b + (up + low) / 2;
        cnt = cnt + 1;
    end
    if cnt == maxIt
        warning('KerMor:svr:Ambiguous_offset',['All coefficients for SVR expansion are bounded.\n'...
            'The offset b cannot be computed sufficiently precise after %d iterations, setting b=%f.'],maxIt,b);
    end
    
else
    % Exact computation of b is possible due to the existance of an alpha_i
    % whose support vector lies on the border of the eps-tube.
    ai = ai(svidx);
    kernelmatrix = kernelmatrix(svidx,svidx);
    % Compute b
    [val,idx] = min(abs(abs(ai)-C/(2*m)));
    % get correct index regarding the complete input set
    index = svidx(idx(1));
    % compute b
    b = fxi(index) - ai' * kernelmatrix(:,idx(1)) + sign(m-index)*eps;
end

end

