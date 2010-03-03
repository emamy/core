function vect = compile_triplevect( model, model_data )

ps = model_data.param_samples;
xi = model_data.snapshots(:,:);

% Repeat times vector as often as different param/input combinations exist
times = repmat(model.times,1,size(model_data.snapshots,3)*size(model_data.snapshots,4));
% Repeat param_samples for each time, then later the result for each input
idx = repmat(reshape(repmat(1:size(ps,2),length(model.times),1),1,[]),1,size(model_data.snapshots,4));
mu = model_data.param_samples(:,idx);
% Compile whole snapshot vector
vect = [times; xi; mu];

end

