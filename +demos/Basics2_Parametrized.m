% Basics2_Parametrized: 
%
% This demo covers:
% - time- and parameter-affine linear system
% - model reduction process for parameterized system with input
% - construction of reduced model
% - simulation of reduced model
% - additional reduction analysis (ModelAnalyzer class)
%
% @author Daniel Wirtz @date 2013-12-17
%
% @new{0,7,dw,2013-12-17} Added this script.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

%% Initialization
% Please refer to demos.Basics1_Linear for elementary linear setup
dim = 1000;
spfac = .2;
% Get a RandStream for reproducable demo
s = RandStream('mt19937ar','Seed',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Simple linear model %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the base model
m = models.BaseFullModel('LinearModelWithParamsDemo');
m.T = .5;
m.dt = 0.01;
m.System.MaxTimestep = m.dt;
% Quick'n dirty. For different solvers see Basics1
m.ODESolver = solvers.ExplEuler;

%% Add some parameters to the system
% This is pretty random! This is for demonstration purposes only.
m.System.addParam('P1 (Diffswitch)',0,'Range',[0 1]);
m.System.addParam('P2 (Exp)',1,'Range',[0 2]);
m.System.addParam('P3 (SomeOtherName)',0,'Range',[-1 1]);

%% Set the basic components
% Define the affine coefficient functions. They can be time- and parameter
% dependent, where the function strings (by convention) have the variables
% `t` and `\mu` available for inline use.
funs = {'t*mu(1)', 'mu(1)*mu(2).^2', 'exp(mu(2))-mu(1)', 'cos(mu(3))+sin(10*mu(1))'};
% Compose an affine-parametric initial value
% Notice: For affine initial values the time parameter will not be used
% (i.e. set to []). Hence, funs{1} would cause an error if used here.
Px0 = dscomponents.AffineInitialValue;
Px0.addMatrix(funs{2},10*(s.rand(dim,1)-.5));
Px0.addMatrix(funs{3},10*(s.rand(dim,1)-.5));
m.System.x0 = Px0;

% Compose an affine-parametric linear core function
% Here we use a diffusion matrix (with neumann boundary) for the first
% component, that will go into the system with `t*\mu(1)`. Hence, the
% diffusion will get stronger over time, and can be controlled via
% parameter one.
Ap = dscomponents.AffLinCoreFun;
e = ones(dim,1);
diff = spdiags([e -2*e e], -1:1, dim, dim);
diff(1,2) = -1; diff(end,end-1) = -1;
Ap.addMatrix(funs{1}, diff);
for k = [2 4]
    Ap.addMatrix(funs{k}, 0.1*(Utils.sprand(dim,dim,spfac,s)-.5));
end
m.System.A = Ap;

% Parameter-dependent output switch between DoFs 1 and 2
Cp = dscomponents.AffLinOutputConv;
Cp.addMatrix('mu(1)',[1; zeros(dim-1,1)]');
Cp.addMatrix('1-mu(1)',[0; 1; zeros(dim-2,1)]');
m.System.C = Cp;

%% Run a simple simulation and plot the results
mu = m.getRandomParam;
[t,y,ct,x] = m.simulate(mu);
fprintf('Simulation time: %g\n',ct);
m.plot(t,y);
m.plotState(t,x);

%% Add inputs
u1 = @(t)-200*exp(-(t-.5).^2/3);
u2 = @(t)200*sin(10*t).*(t>.3);
m.System.Inputs{1} = u1;
m.System.Inputs{2} = u2;
Bp = dscomponents.AffLinInputConv;
Bp.addMatrix('mu(2)/2',[1; zeros(dim-1,1)]);
Bp.addMatrix('1-mu(2)/2',[0; 1; zeros(dim-2,1)]);
m.System.B = Bp;

%% Plot new system
mu = m.getRandomParam;
[~,y] = m.simulate(mu);
[t,y_u1] = m.simulate(mu,1);
m.plot(t,y);
m.plot(t,y_u1);
m.plot(t,y-y_u1);

%% Setup the model reduction
% Choose a parameter sampler
s = sampling.RandomSampler;
s.Seed = 1;
s.Samples = 5;
m.Sampler = s;

% Now we want only the input 1 used as training input (more would result in
% cross-product of param samples and input indices)
m.TrainingInputs = 1;

% Setup the subspace reduction scheme
s = spacereduction.PODGreedy;
% Just go to subspace size 20.
% Otherwise the error tolerance can be set via s.Eps
s.MaxSubspaceSize = 20;
% This causes the trajectory data `x` to be augmented by evaluations
% `A(t,mu)*x`. This can be handy in order to improve the projected A matrix.
% Please see the spacereduction.BaseSpaceReducer class for
% more details and alternatives.
s.IncludeAxData = true;

% This causes the span of the B matrix to be the initial space
% Alternative an initial space can be set using s.InitialSpace
s.IncludeBSpan = true;
m.SpaceReducer = s;

%% Run the offline phase
% All the steps from the offline phase can also be run by using
% m.offlineGenerations.
m.off1_createParamSamples;
m.off2_genTrainingData;
m.off3_computeReducedSpace;
% The following steps are not required in this setting:
% m.off4_genApproximationTrainData;
% m.off5_computeApproximation;
% m.off6_prepareErrorEstimator;
save basics2_param;

%% Simulate & Plot the reduced model
mu = m.getRandomParam;
% This is the MAIN method to create a reduced model from a full model.
r = m.buildReducedModel;
% The models.ReducedModel inherits from BaseModel as well, so that the same
% interface for running simulations can be used.
[~,yr] = r.simulate(mu,1);
[t,y] = m.simulate(mu,1);
m.plot(t,y);
m.plot(t,yr);
m.plot(t,yr-y);

%% Additional analysis tools
ma = ModelAnalyzer(r);
% Gives a short summary of all components of the reduced/full model that
% were involved in the reduction process
ma.plotReductionOverview;
ma.compareRedFull(mu,1);