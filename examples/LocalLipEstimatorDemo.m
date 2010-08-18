function [fm,r,est] = LocalLipEstimatorDemo( dims )
%LOCALLIPESTIMATORDEMO Summary of this function goes here
%   Detailed explanation goes here

%% Global settings
if nargin < 1
    dims = 2000;
end
% Number of support vectors
svNum = 10;
% Strictly positive kernel expansion?
pos_flag = false;
% Uniform expansion or comp-wise separate?
uniform = true;

%% Model settings
fm = models.BaseFullModel;
fm.Verbose = 0;
fm.Name = 'Synthetic kernel based test model 1';

fm.T = 5;
fm.dt = .05;

fm.Approx = [];
fm.Sampler = [];

s = spacereduction.PODReducer;
s.Mode = 'abs';
s.Value = 1;
fm.SpaceReducer = s;

%this.ODESolver = solvers.MLWrapper(@ode45);
fm.ODESolver = solvers.ExplEuler;

%% Core function
cf = dscomponents.CompwiseKernelCoreFun;
cf.SystemKernel = kernels.GaussKernel(dims*100);
cf.TimeKernel = kernels.NoKernel;
cf.ParamKernel = kernels.NoKernel;
cf.snData.xi = repmat(linspace(-20,20,svNum),dims,1);
cf.snData.ti = [];
cf.snData.mui = [];

% Function coefficients
offset = .5;
if pos_flag
    offset = 0;
end
% Create coefficients
if uniform
    ai = (rand(1,svNum)-offset);
    cf.Ma = repmat(ai,dims,1);
else
    cf.Ma = (rand(dims,svNum)-offset);
end

%% System settings
s = models.BaseDynSystem;

%s.x0 = @(mu)ones(dims,1)*.5;
s.x0 = @(mu)rand(dims,1)*.5;

s.f = cf;

%s.Inputs{1} = @(t)0;

fm.System = s;

%% Generation
fm.offlineGenerations;
r = fm.buildReducedModel;

%% Error estimators
est = struct;
est(1).Name = 'Full error';
est(1).Estimator = error.DefaultEstimator(r);
est(1).Estimator.Enabled = true;
est(2).Name = 'Global Lipschitz estimator';
est(2).Estimator = error.GlobalLipKernelEstimator(r);
est(3).Name = 'Local Lipschitz estimator: getLocalGradientLipschitz';
est(3).Estimator = error.LocalLipKernelEstimator(r);
k = cf.SystemKernel;
est(3).Estimator.KernelLipschitzFcn = @k.getLocalGradientLipschitz;
est(4).Name = 'Local Lipschitz estimator: getLocalSecantLipschitz';
est(4).Estimator = error.LocalLipKernelEstimator(r);
est(4).Estimator.KernelLipschitzFcn = @k.getLocalSecantLipschitz;
est(5).Name = 'Local Lipschitz estimator: getImprovedLocalSecantLipschitz';
est(5).Estimator = error.LocalLipKernelEstimator(r);
est(5).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;

%% Simulations
num = 5;
times = zeros(num,1);
errs = zeros(num,length(fm.Times));
figure;
hold on;
for idx = 1:num
    r.ErrorEstimator = est(idx).Estimator;
    [t,y,times(idx)] = r.simulate;
    errs(idx,:) = r.ErrorEstimator.LastError;
end
%T = repmat(t,num,1) ./ repmat(times,1,length(fm.Times));
%% Plot
figure(1);
a = cell(1,5);
[a{:}] = est(:).Name;

subplot(1,2,1);
plot(fm.Times,errs);
xlabel('Time');
ylabel('Error');
legend(a);

subplot(1,2,2);
style = 's';
plot(errs(1,end),times(1),style,...
    errs(2,end),times(2),style,...
    errs(3,end),times(3),style,...
    errs(4,end),times(4),style,...
    errs(5,end),times(5),style,...
    'MarkerSize',7);
%txt = [repmat('(',num,1) num2str(times) repmat(', ',num,1) num2str(errs(:,end)) repmat('s)',num,1)];
%text(times,errs(:,end),txt);
xlabel('e(T)');
ylabel('Comp. time');
legend(a);
%pause;
%close(h);

end

