function model_data = gen_snapshots( model, model_data )
%GENERATE_SNAPSHOTS Summary of this function goes here
%   Detailed explanation goes here

ps = model_data.param_samples;

% Need minimum "one" for loops.
num_samples = max(1,size(ps,2));
num_inputs = max(1,length(model.system.inputs));
num_times = length(model.times);

% Compute system dimension using x0.
mu = [];
if ~isempty(model.system.params)
    mu = zeros(length(model.system.params),1);
end
dims = length(model.system.x0(mu));
% Initialize snapshot array
snapshots = zeros(dims, length(model.times), num_samples, num_inputs);
f_val = zeros(dims, length(model.times), num_samples, num_inputs);

try
    wh = waitbar(0,'Initializing snapshot generation...');
    cnt = 1;
    % Iterate through all input functions
    for inidx = 1:num_inputs
        % Iterate through all parameter samples
        for pidx = 1:num_samples
            
            % Check for no parameters
            if isempty(ps)
                mu = [];
            else
                mu = ps(:,pidx);
            end
            % Check for no inputs
            if isempty(model.system.inputs)
                inputidx = [];
            else
                inputidx = inidx;
            end
            
            % Display
            perc = cnt/(num_inputs*num_samples);
            waitbar(perc,wh,sprintf('Generating snapshots ... %2.0f %%',perc*100));
            cnt=cnt+1;
            
            % Get ODE function (general function)
            [t,x] = compute_trajectory(model, mu, inputidx);
            % Store each solution/trajectory in rows=dimension, column=timestep
            x = x';
            % Assign snapshot value
            snapshots(:,:,pidx,inidx) = x;
            % Evaluate f at those points
            for tidx=1:num_times
                fx = model.system.f(x(:,tidx),model.times(tidx),mu);
                f_val(:,tidx,pidx,inidx) = fx;
            end
        end
    end
    
catch ME
    close(wh);
    rethrow(ME);
end

model_data.snapshots = snapshots;
model_data.f_values = f_val;

if ishandle(wh)
    close(wh); 
end
