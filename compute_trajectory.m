function [ t,x ] = compute_trajectory(model, mu, inidx)
%COMPUTE_TRAJECTORY Computes a trajectory of a model for specific mu and input
%
% Params:
%
%

odefun = gen_odefun(model.system, mu, inidx);
x0 = model.system.x0(mu);

% Solve ODE
opts = [];
if ~isinf(model.system.max_timestep)
    opts = odeset('MaxStep',model.system.max_timestep, 'InitialStep',.5*model.system.max_timestep);
end

[t,x] = model.odesolver(odefun, model.times, x0, opts);
end