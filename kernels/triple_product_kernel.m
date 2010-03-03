function matrix = triple_product_kernel(tri1,tri2, xdim, timekernel, spacekernel, paramkernel)
    % Arguments are matrices with t,x and mu samples in each column
    t1 = tri1(1,:);
    x1 = tri1(2:xdim+1,:);
    mu1 = tri1(xdim+2:end,:);
    t2 = tri2(1,:);
    x2 = tri2(2:xdim+1,:);
    mu2 = tri2(xdim+2:end,:);

    matrix = timekernel(t1,t2) .* spacekernel(x1,x2) .* paramkernel(mu1,mu2);
end

