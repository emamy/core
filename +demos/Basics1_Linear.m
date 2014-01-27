% Basics_Linear: Demo file for KerMor linear models and reduction. 
%
% Covers:
% - Linear system with inputs and outputs
% - Basic plotting of results
% - Different solver selection
% - Mass matrices
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

% Metadata
dim = 4;
s = RandStream('mt19937ar','Seed',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Simple linear model %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the base model
m = models.BaseFullModel('LinearModelDemo');
% Max simulation time
m.T = 1.5;
% The desired timestepping (output)
m.dt = 0.01;
% The internal maximum timestep for any solver (details later)
m.System.MaxTimestep = m.dt;
disp(m);

%% Set the basic components
% Only a core linear matrix A and an initial condition are needed here.
m.System.A = dscomponents.LinearCoreFun(s.rand(dim,dim));
m.System.x0 = dscomponents.ConstInitialValue(10*(s.rand(dim,1)-.5));

%% Run a simple simulation and plot the results
[t,y,ct,x] = m.simulate;
fprintf('Simulation time: %g\n',ct);
% Plot the system output
m.plot(t,y);
% Plot the state space - here the same, as no specific output conversion is
% set
m.plotState(t,x);

%% Modify the inputs and outputs
% Now we only want the first DoF as output
m.System.C = dscomponents.LinearOutputConv([1; zeros(dim-1,1)]');
% And we also add some system inputs!
% This is split into B*u(t), so we need both components.
u1 = @(t)-200*exp(-(t-.5).^2/3);
u2 = @(t)200*sin(10*t).*(t>.3);
m.System.Inputs{1} = u1;
m.System.Inputs{2} = u2;
figure;
plot(m.Times,u1(m.Times),'r',m.Times,u2(m.Times),'b');
legend('u1','u2');
% Lets map the inputs to the second DoF
m.System.B = dscomponents.LinearInputConv([0; 1; zeros(dim-2,1)]);

%% Run simulation and plot the results
% The first argument is a parameter by convention - so leave [] for input
% selection.
[t,y,~,x] = m.simulate([],1);
m.plot(t,y);
m.plotState(t,x);

[t,y,~,x] = m.simulate([],2);
m.plot(t,y);
m.plotState(t,x);

%% Lets use a different solver
% Choose explicit Heun's method
m.ODESolver = solvers.Heun;
[t,y_heun] = m.simulate([],2);
m.plot(t,y_heun);

% Choose the MatLab-builtin ode15i implicit solver
m.ODESolver = solvers.MLode15i;
[t,y_impl] = m.simulate([],2);
m.plot(t,y_impl);
figure;
semilogy(t,abs(y_impl-y_heun));

%% Adding a mass matrix E
% Choose a positive definite matrix
Mmat = diag(s.randi(20,1,dim)) + s.rand(dim,dim);
m.System.M = dscomponents.ConstMassMatrix(Mmat);
m.ODESolver = solvers.Heun;
[t,y] = m.simulate([],1);
m.plot(t,y);

m.ODESolver = solvers.MLode15i;
[t,y2] = m.simulate([],1);
m.plot(t,y2);