% Basics5_KernelsAndApprox: 
%
% This covers:
% - Basic concepts regarding kernel expansions
% - VKOGA approximation
% - Visualization
%
% @author Daniel Wirtz @date 2013-12-18
%
% @new{0,7,dw,2013-12-18} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

%% Create a basic kernel expansion manually
lw = 2;
% Centers
c = [-1 .5 1.2];
% Coefficients
a = [.5 -.2 .8];
% Get a Gaussian kernel with gamma=1
k = kernels.GaussKernel(1);
% Get the overall linear combination class and set properties
f = kernels.KernelExpansion;
f.Kernel = k;
f.Centers.xi = c;
f.Ma = a;

%% Plot the expansion and its centers/single elements
% More plot functions illustrating kernels can be found at
% kernels.Wendland.test_WendlandKernel and
% kernels.GaussKernel.test_GaussKernel
pm = PlotManager;
h = pm.nextPlot('rkhs_demo','Example kernel expansion','x','f(x)');
x = -3:.01:3;
plot(h,x,f.evaluate(x),'r','LineWidth',lw);
hold(h,'on');
plot(h,x,0,'k','LineWidth',lw);
for i = 1:length(c)
    f.Centers.xi = c(i);
    f.Ma = a(i);
    plot(h,x,f.evaluate(x),'r--');
    plot(h,c(i),0,'b.','MarkerSize',20);
    plot(h,[c(i) c(i)+eps],[0 a(i)],'k--');
end
pm.done;

%% Run a VKOGA algorithm
load +demos\spine_training_data_set.mat;

% Get VKOGA instance
alg = approx.algorithms.VKOGA;
% Set maximal size
alg.MaxExpansionSize = 600;
% Choose variant
alg.UsefPGreedy = false;

kexp = kernels.KernelExpansion;
nG = 10;
k = [1 3];
kexp.Kernel = kernels.Wendland;

% Create configuration for VKOGA
gammas = linspace(15,70,nG);
comb = Utils.createCombinations(gammas,k);
wc = kernels.config.WendlandConfig('G',comb(1,:),'S',comb(2,:));
wc.Dimension = 3;
ec = kernels.config.ExpansionConfig;
ec.StateConfig = wc;
alg.ExpConfig = ec;

%% Execute computations
alg.computeApproximation(kexp, atd);
kexp_dbase = kexp.toTranslateBase;
save basics5_kernelsandapprox;
