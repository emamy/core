function model_data = gen_reduced_space( model, model_data )
%GEN_REDUCED_SPACE Summary of this function goes here
%   Detailed explanation goes here

% Get snapshots as col-vect x instance matrix
if ~isfield(model_data,'snapshots')
    error('Field ''snapshots'' missing in model_data. Called gen_snapshots(model, param_samples) yet?');
end

switch model.reduction.mode
    case 'POD'
        % Call POD reducer with POD settings
        model_data = reduce_pod(model, model_data);
end


end

