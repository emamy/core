% Basics3_Nonlinear: 
%
% - A nonlinar (burgers) model
% - nonlinear MOR with DEIM-nonlinearity approximation
% - 
%
% @author Daniel Wirtz @date 2013-12-18
%
% @new{0,7,dw,2013-12-18} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

% Please see Basics1 and Basics2 for more elementary setups. For the
% burgers model, most of the setup has been moved into class constructors
% of classes inheriting from models.BaseFullModel and models.BaseFirstOrderSystem.

dim = 500;

%% Model and system setup
% The original burgers model has a nonzero initial condition and no inputs.
% we add the same things here as also been used in the a-posteriori error
% estimation paper.
m = models.burgers.Burgers(dim,2);
% This solver solves the linear part implicitly and the nonlinear party
% explicitly.
m.ODESolver = solvers.SemiImplicitEuler(m);
% This determines the plot angle for the burgers plots.
m.PlotAzEl = [-49 34];
m.SaveTag = sprintf('burgers_d%d_fx1_bs1',m.Dimension);

s = m.System;
% Zero initial value
s.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
% Add Inputs
x = linspace(m.Omega(1),m.Omega(2),dim+2);
x = x(2:end-1);
pos1 = logical((x >= .1) .* (x <= .3));
pos2 = logical((x >= .6) .* (x <= .7));
s.Inputs{1} = @(t)[sin(2*t*pi); (t>.2)*(t<.4)];
B = zeros(dim,2);
B(pos1,1) = 4*exp(-((x(pos1)-.2)/.03).^2);
B(pos2,2) = 4;
s.B = dscomponents.LinearInputConv(B);

%% Setup reduction
m.TrainingInputs = 1;

% Sampling - log-grid
p = m.System.Params(1);
p.Range = [0.01, 0.06];
p.Desired = 20;
p.Spacing = 'log';
m.Sampler = sampling.GridSampler;

% Space reduction
% m.ComputeTrajectoryFxiData = true;
p = spacereduction.PODReducer;
% p.IncludeTrajectoryFxiData = true;
p.IncludeFiniteDifferences = false;
p.IncludeBSpan = true;
p.Mode = 'abs';
p.Value = 50;
m.SpaceReducer = p;
            
% Approximation of nonlinearity: Choose DEIM method here
a = approx.DEIM(m.System);
a.MaxOrder = 100;
m.Approx = a;

% Run offline phase
m.offlineGenerations;
save basic3_nonlinear;
%load basic3_nonlinear;

%% Build reduced model and do some analysis
% Details on construction of reduced models see Basics2!
r = m.buildReducedModel;
% Set the reduced DEIM approximation order to 25, and 10 orders for the
% error estimator.
r.System.f.Order = [25 10];
mu = m.getRandomParam(1,1);
ma = ModelAnalyzer(r);
ma.compareRedFull(mu,1);
ma.analyzeError(mu,1);

%% Goodie: start the DEIM error estimator analyzer!
DEIMEstimatorAnalyzer(r);