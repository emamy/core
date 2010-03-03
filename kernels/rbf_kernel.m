function K = rbf_kernel( X,Y,gamma )
%RBF_KERNEL Summary of this function goes here
%   Detailed explanation goes here

n1sq = sum(X.^2,1);
n1 = size(X,2);

if isempty(Y);
    n2sq = n1sq;
    n2 = n1;
else
    n2sq = sum(Y.^2,1);
    n2 = size(Y,2);
end;
K = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*X'*Y;
K(K<0) = 0;
K = exp(-gamma*K);
end