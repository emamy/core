function find_vec_in_matrix_speedtest
% Created 2011-04-01 for test purposes of finding a vector in a matrix.

dim1 = 500;
dim2 = 6000;
A = rand(dim1,dim2);
b = A(:,round(dim2/2));
b(1) = 200;

time = tic;
%f = @(x,y)(x-y).^2;
hlp = bsxfun(@eq,A,b);
ind = find(sum(hlp,1) == dim1,1)
time1 = toc(time);

time = tic;
ind = find(sum(A == repmat(b,1,dim2)) == dim1,1)
time2 = toc(time);

time = tic;
ind = strfind(reshape(A,1,[]),b');
ind = round((ind+dim1-1)/dim1)
time3 = toc(time);

fprintf('Dims: %dx%d, bsxfun time: %f, repmat time: %f, strfind time: %f\n',dim1,dim2,time1,time2,time3);

end