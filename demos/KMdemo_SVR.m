function out = KMdemo_SVR(in)
% KMdemo_SVR: Demonstrates the ScalarEpsSVR_SMO class
%
%
%
% @author Daniel Wirtz @date 2011-10-05
%
% @new{0,5,dw,2011-10-05} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

% Performs a test of this class

if nargin == 0
    version = 1;
end

x = -5:.1:5;
fx = sinc(x)+.2*x;
%             fx = sinc(x);

svr = general.regression.ScalarEpsSVR_SMO;
svr.Version = version;
%svr.Eps = 0.073648;
svr.Eps = .1;
svr.Lambda = 1/20;%1/20; % i.e. C=10 as in ScalarEpsSVR
svr.Vis = 1;

%kernel = kernels.PolyKernel(7);
%kernel = kernels.LinearKernel;
kernel = kernels.GaussKernel(.8);
svr.K = data.MemoryKernelMatrix(kernel.evaluate(x,x));

[ai, svidx] = svr.computeKernelCoefficients(fx,[]);
sv = x(:,svidx);
svfun = @(x)ai'*(kernel.evaluate(x,sv)');

fsvr = svfun(x);

fdiff = abs(fsvr(svidx)-fx(svidx));
errors = find(fdiff > 1.01*svr.Eps);
res = isempty(errors);

% Plot approximated function
figure;
plot(x,fx,'r',x,[fx-svr.Eps; fx+svr.Eps],'r--');
hold on;
plot(x,fsvr,'b',x,[fsvr-svr.Eps; fsvr+svr.Eps],'b--');
skipped = setdiff(1:length(x),svidx);
plot(sv,fx(svidx),'.r','MarkerSize',20);
plot(x(skipped),fx(skipped),'xr');

if ~res
    plot(x(svidx(errors)),fx(svidx(errors)),'blackx','LineWidth',4);
end

tit = sprintf('#SV=%d, eps=%f',length(svidx),svr.Eps);
title(tit);
disp(tit);
hold off;

end