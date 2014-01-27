function rsPenaltyFunTest
% rsPenaltyFunTest:
%
% Checks if the penlized newton function to find `r_s` always finds the correct root and the newton
% penalization polynomials dont introduce any root above `r_m` for bad bell function configuration.
%
% @author Daniel Wirtz @date 2011-05-30
%
% @new{0,4,dw,2011-05-30} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

k = kernels.GaussKernel(1);
figure(1);

g = logspace(-3,2,30);
for gamma = g
    k.Gamma = gamma;
    
    x = linspace(0,2*k.rm,200);
    
    fx = k.evaluateScalar(x);
    
    nx = k.evaluateD1(x) - (k.evaluateScalar(x)-1)./x;
    nx2 = k.evaluateD1(x) - (k.evaluateScalar(x)-k.evaluateScalar(k.r0/2))./(x-k.r0/2);
    nx3 = k.evaluateD1(x) - (k.evaluateScalar(x)-k.evaluateScalar(k.r0))./(x-k.r0);
    
    plot(x,fx,'r',x,nx,'b',x,nx2,'g',x,[nx3; 3*[0 diff(nx3)]],'m',[k.rm k.rm],[0 k.evaluateD1(k.rm) - (k.evaluateScalar(k.rm)-1)/k.rm],'black',x,0,'black');
    axis tight;
    pause;
end