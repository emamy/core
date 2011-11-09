function times = KMDemo_AffParamMatrixSpeed(n)
% KMDemo_AffParamMatrixSpeed: 
%
%
%
% @author Daniel Wirtz @date 2011-10-25
%
% @new{0,5,dw,2011-10-25} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

if nargin < 1
    n = 5000;
end

d = 50;

str = {'4*t','4*t + mu(1)','exp(sum(mu))','sin(t)','exp(-t)*mu(2)','t*mu(2)','cos(t)','exp(-mu(2))'};
theta{1} = @(t,mu)4*t;
theta{2} = @(t,mu)4*t + mu(1);
theta{3} = @(t,mu)exp(sum(mu));
theta{4} = @(t,mu)sin(t);
theta{5} = @(t,mu)exp(-t)*mu(2);
theta{6} = @(t,mu)t*mu(2);
theta{7} = @(t,mu)cos(t);
theta{8} = @(t,mu)exp(-mu(2));
nc = length(theta);

am = general.AffParamMatrix;

for i=1:2
    mat{i} = rand(d,d);%#ok
    am.addMatrix(str{i}, mat{i});
end

v = rand(d,1);
params = rand(3,n);

%% Linear
% t = tic;
% for i=1:n
%     coeff = zeros(1,5);
%     for k=1:nc
%         fun = theta{k};
%         coeff(k) = fun(params(1,i),params(2:3,i));
%     end
%     dummy = lincomb_sequence(mat,coeff)*v;%#ok
% end
% t1 = toc(t);
% 
% t = tic;
% for i=1:n
%     dummy = am.compose(params(1,i),params(2:3,i))*v;%#ok
% end
% t2 = toc(t);
% 
% times = [t1 t2];
times = [];

%% Quadratic
t = tic;
for i = 1:nc
    for j = 1:nc
        mat2{i,j} = mat{i}*mat{j};%#ok
    end
end
for i=1:n
    coeff = zeros(1,5);
    for k=1:nc
        fun = theta{k};
        coeff(k) = fun(params(1,i),params(2:3,i));
    end
    dummy = lincomb_sequence2(mat2,coeff,coeff)*v;%#ok
end
t1 = toc(t);

t = tic;
am = am*am;
for i=1:n
    dummy = am.compose(params(1,i),params(2:3,i))*v;%#ok
    %dummy = am(params(1,i),params(2:3,i))*v;%#ok
end
t2 = toc(t);

times = [times t1 t2];
end