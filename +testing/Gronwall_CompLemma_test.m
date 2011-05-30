function Gronwall_CompLemma_test
% Gronwall_CompLemma_test: 
%
% Tests the estimation sharpness of the Gronwall's Lemma and the Comparison Lemma
%
% @author Daniel Wirtz @date 2011-05-25
%
% @new{0,4,dw,2011-05-25} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

T = 5;
t = 0:.01:T;

e0 = 1;
aG = @(t)2*t + e0;
aC = @(t)2; % = aG'
b = @(t)2*(T-t)/T;

odefunGW = @(t,y)b(t)*(aG(t) + y);
odefunCL = @(t,y)b(t)*y + aC(t);

opts = odeset(odeset,'RelTol',1e-8,'AbsTol',1e-9);
[t,yG] = ode113(odefunGW,t,0,opts);
[t,uC] = ode113(odefunCL,t,e0,opts);

uG = aG(t) + yG;

figure(1);
subplot(1,2,1);
semilogy(t,uG,'r',t,uC,'b');
legend('\Delta_G(t)','\Delta_C(t)','Location','NorthWest');
subplot(1,2,2);
%plot(t,(uC-uG),'r');
semilogy(t,max(0,uC-uG),'r',t,max(0,-(uC-uG)),'b');
legend('\Delta_G(t)-\Delta_C(t)');
fprintf('Factor G/L=%f\n',uG(end)/uC(end));
end