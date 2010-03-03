function testsvr
%TESTSVR Summary of this function goes here
%   Detailed explanation goes here

x = -2:.1:5;
fx = sinc(x);

figure(1);
eps = .1;
plot(x,fx,'r',x,[fx-eps; fx+eps],'r--');

kernel = get_rbf_kernelfun(2);

[ai,b,svidx] = scalar_svr(fx,kernel(x,x),eps,100000);
sv = x(svidx);
svfun = @(x)ai'*kernel(sv,x) + b;

fsvr = svfun(x);

hold on;

% Plot approximated function
plot(x,fsvr,'b',x,[fsvr-eps; fsvr+eps],'b--');
skipped = setdiff(1:length(x),svidx);
plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');

hold off;

end

