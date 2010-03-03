function kernelfun = get_rbf_kernelfun( gamma )
%GET_RBF_KERNEL Summary of this function goes here
%   Detailed explanation goes here

kernelfun = @(X,Y)rbf_kernel(X,Y,gamma);

end

