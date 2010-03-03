function K = poly_kernel( X,Y,deg )
%POLY_KERNEL Summary of this function goes here
%   Detailed explanation goes here

K = X'*Y.^deg;

end

