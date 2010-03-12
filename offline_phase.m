function [model_data, reduced_data] = offline_phase( model )
%OFFLINE_PHASE Summary of this function goes here
%   Detailed explanation goes here

% Create detailed model data
model_data = gen_model_data(model);
% Compute reduced data
reduced_data = gen_reduced_data(model, model_data);

end

