function s = testsettings
%TESTSETTINGS KerMor test settings collection
%
%   Common structure from which unit tests throughout the framework should
%   use settings and functions.
%   This is intended to ensure tests use the same data base to get some
%   kind of comparability.
%
% @author Daniel Wirtz @date 16.03.2010
s = struct;

%% Test sizes
% What dimension?
s.testdim = 10;

% How many params? (to use from the ones defined below)
s.testparams = 2;
s.testinputs = 1;

%% Model settings
s.m = models.BaseFullModel;
s.m.T = 3;
s.m.dt = .2;
s.m.Sampler = sampling.RandomSampler;
s.m.Sampler.Samples = 10;
s.m.System.MaxTimestep = s.m.dt;
s.m.ODESolver = solvers.ExplEuler;

a = approx.KernelApprox;
ap = approx.algorithms.VKOGA;
ec = kernels.config.ParamTimeExpansionConfig;
ec.StateConfig = kernels.config.GaussConfig('G',1:2);
ec.ParamConfig = kernels.config.GaussConfig('G',.1:.1:.3);
ap.ExpConfig = ec;
ap.MaxExpansionSize = 5; % Keep test run short!
a.Algorithm = ap;
%a.CoeffComp = general.regression.ScalarEpsSVR;
s.TimeKernel = kernels.NoKernel;
%a.TimeKernel = kernels.GaussKernel(2);
s.Kernel = kernels.GaussKernel(2);
s.ParamKernel = kernels.GaussKernel(2);
s.m.Approx = a;

%% Dynamical System settings/functions
s.Inputs{1} = @(t)1; % Function 1: Constant 1
s.Inputs{2} = @(t)sin(4*t); % Function 2: some sin(t)
s.Inputs = s.Inputs(1:s.testinputs);

% Used Parameter Space
s.params(1) = struct('Name', 'P1', 'MinVal', -1, 'MaxVal', 1, 'Desired', 10);
s.params(2) = struct('Name', 'P2', 'MinVal', 2, 'MaxVal', 3, 'Desired', 5);
s.params(3) = struct('Name', 'P3', 'MinVal', 0, 'MaxVal', 10, 'Desired', 6);
s.params = s.params(1:s.testparams);

% Input conversion
muidx = randi(s.testparams,s.testdim,1);
s.B = @(t,mu)ones(s.testdim,1);
s.B_p = @(t,mu)ones(s.testdim,1).*mu(muidx);

% Initials
musel = (rand(s.testdim,1)<.5);
s.x0 = dscomponents.ConstInitialValue(rand(s.testdim,1));
av = dscomponents.AffineInitialValue;
av.addMatrix(['mu(' num2str(randi(s.testparams)) ')'],musel .* rand(s.testdim,1));
av.addMatrix('.5 + mu(1)',5*rand(s.testdim,1));
s.x0_p = av;

% Linear functions
muidx = randi(s.testparams,s.testdim,1);
A = rand(s.testdim,s.testdim);
s.flin = @(x,t,mu)A*x;
s.flin_p = @(x,t,mu)(A + diag(mu(muidx)))*x;

% Nonlinear functions
randx = randi(s.testdim);
randmu = randi(s.testparams);
s.fnlin_p = @(x,t,mu)(.5+t/2)*sin(x) + x(randx)*mu(randmu);
s.fnlin = @(x,t,mu)(.5+t/2).*sin(x);

end

