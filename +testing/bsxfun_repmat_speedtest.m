function bsxfun_repmat_speedtest
% Created 2011-04-01 for test purposes of finding all distances from one
% vector to a list of vectors
%
% Result: repmat is faster!

dim1 = 50;
dim2 = 80;
A = rand(dim1,dim2);
b = A(:,round(dim2/2));
%b(1) = 200;

time = tic;
hlp = bsxfun(@(x,y)(x-y).^2,A,b);
dist = sqrt(sum(hlp,1));
time1 = toc(time);

time = tic;
B = repmat(b,1,dim2);
dist = sqrt(sum((A - B).^2,1));
time2 = toc(time);

fprintf('Dims: %dx%d, bsxfun time: %f, repmat time: %f\n',dim1,dim2,time1,time2);

end