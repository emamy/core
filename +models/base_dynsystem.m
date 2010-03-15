function dynsys = base_dynsystem
%base_dynsystem KerMor base dynamical system struct
%
% Use this as starting point for specific dynamical systems.

% Overall dynamical system struct
dynsys = struct;

%% General Settings
% This setting is important if the f function is in fact a
% space-discretized PDE where timestep constraints (CFL) apply, for
% example. 
% Default: Inf (any timestep possible)
dynsys.max_timestep = Inf;
% This struct shall contain all system-specific information which is
% irrelevant for the KerMor framework.
dynsys.specific = struct;

%% Initial value function
% Default: Zero initial value, dependent on mu
% Format: Function handle that takes a parameter mu and returns a col-vector
dynsys.x0 = @(mu)0;

%% System Inputs
% Default: No input conversion
% Format: Function handle, dependent on t,mu
dynsys.B = @(t,mu)0;
%dynsys.B = struct('coeff',@(t,mu)0,'B',ones(2,2));
% Default: No inputs. Has to be a cell of function handles
% Format: Function handle each, dependent on t
dynsys.inputs = [];

%% System parameter Space
% Default: No parameters.
% Format: params is a struct array with each struct in the format of
%         struct('Name',<name>,'MinVal',<min>,'MaxVal',<max>,'Desired',<num>)
dynsys.params = struct('Name',{},'MinVal',{},'MaxVal',{},'Desired',{});

%% Nonlinearity definitions
% Default: Zero
% Type: function handle, dependent on x,t,mu
dynsys.f = @(x,t,mu)0;

%% System outputs
% Default: Empty
% Type: Function handle, dependent on t,mu
dynsys.C = @(t,mu)0;

end