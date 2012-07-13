function t = BinaryvsMatSaveSpeed
% BinaryvsMatSaveSpeed: 
%
%
%
% @author Daniel Wirtz @date 2012-07-09
%
% @new{0,6,dw,2012-07-09} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

runs = 2;
t = zeros(3,runs);
for k = 1:runs
    %A = rand(3355,2000);
    A = full(sprand(5000,10000,.4));
    ti = tic;
    save A.mat A;
    t(1,k) = toc(ti);
    
    ti = tic;
    save Av4.mat A -v4;
    t(2,k) = toc(ti);
    
    ti = tic;
    save Av6.mat A -v6;
    t(2,k) = toc(ti);
    
%     ti = tic;
%     f = fopen('A.bin','w+');
%     fwrite(f,A,'double');
%     fclose(f);
%     t(3,k) = toc(ti);
end
t = mean(t,2);
end