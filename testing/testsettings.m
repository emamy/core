function s = testsettings
%TESTSETTINGS Summary of this function goes here
%   Detailed explanation goes here

s = struct;

%% Test sizes
% What dimension?
s.testdim = 5;
% How many params? (to use from the ones defined below)
s.testparams = 1;
s.testinputs = 0;

%% Model settings
s.times = 0:.1:1;
s.mode = 'rand';
s.samples = 10;


%% Dynamical System settings/functions

s.inputs{1} = @(t)1; % Function 1: Constant 1
s.inputs{2} = @(t)sin(4*t); % Function 2: some sin(t)
s.inputs = s.inputs(1:s.testinputs);

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
muidx = randi(s.testparams,s.testdim,1);
musel = (rand(s.testdim,1)<.5);
s.x0 = @(mu)rand(s.testdim,1);
s.x0_p = @(mu)rand(s.testdim,1)+musel.*mu(muidx);

% Linear functions
muidx = randi(s.testparams,s.testdim,1);
A = rand(s.testdim,s.testdim);
s.flin = @(x,t,mu)A*x;
s.flin_p = @(x,t,mu)(A + diag(mu(muidx)))*x;

% Nonlinear functions
randx = randi(s.testdim);
randmu = randi(s.testparams);
s.fnlin_p = @(x,t,mu)(.5+t)*sin(x) + x(randx)*mu(randmu);
s.fnlin = @(x,t,mu)(.5+t)*sin(x);

end

