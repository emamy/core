function dynsys = test_dynsys
%TEST_DYNSYS Dynamical System for testing purposes

ts = testsettings;

%% Load basic dynamical system
dynsys = base_dynsystem;

%% Initial value function
dynsys.x0 = ts.x0_p;

%% System Inputs
dynsys.B = ts.B_p;
dynsys.inputs = ts.inputs;

%% System parameter Space
dynsys.params = ts.params;

%% Nonlinearity definitions
dynsys.f = ts.fnlin_p;

end

