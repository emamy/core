function LoopTest(d)
% LoopTest: A quick test to show how slow looping in MatLab can be.
%
% @author Daniel Wirtz @date 2011-05-09
%
% @new{0,4,dw,2011-05-09} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

if nargin < 1
    d = 200;
end
n = 50;
idx = 1:d;
A = zeros(size(d,d));
fprintf('Starting loop speed-test with d=%d and n=%d runs...\n',d,n);
times = zeros(4,n);
for r = 1:n
    t = tic;
    for i = idx
        for j = idx
            A(i,j) = exp(-i+j*.3) - sin(i*j*pi);
        end
    end
    times(1,r) = toc(t);

    t = tic;
    [x,y] = meshgrid(idx,idx);
    all = [reshape(x,1,[]); reshape(y,1,[])];
    for i = 1:d*d
        A(all(1,i),all(2,i)) = exp(-all(1,i)+all(2,i)*.3) - sin(all(1,i)*all(2,i)*pi);
    end
    times(2,r) = toc(t);

    t = tic;
    [x,y] = meshgrid(idx,idx);
    all = [reshape(x,1,[]); reshape(y,1,[])];
    for k = 1:d*d
        i = all(1,k);
        j = all(2,k);
        A(i,j) = exp(-i + j*.3) - sin(i * j * pi);
    end
    times(3,r) = toc(t);
    
    t = tic;
    [y,x] = meshgrid(idx,idx);
    A = exp(-x + y*.3) - sin(x .* y * pi);
    times(4,r) = toc(t);
end
t = sum(times,2) / n;
disp(times);
fprintf('%fs\n',t);
end