function model_data = offline_phase( model )
%OFFLINE_PHASE Summary of this function goes here
%   Detailed explanation goes here

% Sampling of parameter space (if set)
model_data = gen_param_samples(model);
% Compute system snapshots / f_i's
% (adds model_data.snapshots, model_data.f_values)
model_data = gen_snapshots(model, model_data);
% Compute reduced space
% (adds model_data.V)
model_data = gen_reduced_space(model, model_data);
% Compute approximation of f
% (adds model_data.suppvect, model_data.approx_data)
model_data = gen_approximation(model, model_data);

end

