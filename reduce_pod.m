function model_data = reduce_pod(model, model_data)
%REDUCE_POD Summary of this function goes here
%   Detailed explanation goes here

%% Validation
podsettings = model.reduction.POD;
if ~isfield(podsettings,'mode')
    error('''mode'' field in model.system.POD is missing.');
elseif ~isfield(podsettings,'value')
    error('''value'' field in model.system.POD is missing.');
elseif ~any(strcmpi(podsettings.mode,{'sign','eps','rel','abs'}))
    error(['unknown POD reduction mode: ''' podsettings.mode '''']);
end

%% Preparation
% Collapse parameter samples and input dimension
sn = model_data.snapshots(:,:);

A = sn'*sn;

useeig = false;
%% Dimension checks
if strcmpi(podsettings.mode,'rel')
    target_dim = round(size(sn,1)*podsettings.value);
    if target_dim > size(A,2)
        warning('KerMor:POD:mode_rel','Reduced space cant be bigger (%d) than sample size (%d). Using sample size.',target_dim,size(A,2));
        useeig = true;
    end
elseif strcmpi(podsettings.mode,'abs')
    target_dim = podsettings.value;
    if target_dim > size(A,2)
        warning('KerMor:POD:mode_abs','''model.reduction.POD.value'' cant be bigger (%d) than sample size (%d). Setting to sample size.',podsettings.value,size(A,2));
        podsettings.value = size(A,2);
        useeig = true;
    end
end
%% Actual computation
if any(strcmpi(podsettings.mode,{'sign','eps'})) || useeig
    % compute full eigenvalues for sign/eps mode or if other modes result
    % in full space, too. (efficiency)
    [ev,ew] = eig(A);
    % Invert
    ev = flipdim(ev,2);
    ew = flipdim(diag(ew),1);
else
    % Forward verbose setting to eigs fcn
    opts.disp = model.info.verbose;
    [ev, ew] = eigs(A,target_dim,'lm',opts);
    ew = diag(ew);
end
%% Reduction for modes 'sign' and 'eps'
if strcmpi(podsettings.mode,'sign')
    sig = ew >= ew(1)*podsettings.value;
    % Reduce
    ev = ev(:,sig);
    ew = ew(sig);
elseif strcmpi(podsettings.mode,'eps')
    sig = ew >= podsettings.value;
    % Reduce
    ev = ev(:,sig);
    ew = ew(sig);
end

% switch lower(podsettings.mode)
%     % Significance choice: Only use eigenvectors whose eigenvalues are
%     % greater than 'value' percent of the largest eigenvalue
%     case 'sign'
%         [ev,ew] = eig(A);
%         % Invert
%         ev = flipdim(ev,2);
%         ew = flipdim(diag(ew),1);
%         sig = ew >= ew(1)*podsettings.value;
%         % Reduce
%         ev = ev(:,sig);
%         ew = ew(sig);
%     % Use evectors with evalues larger than 'value'
%     case 'eps'
%         [ev,ew] = eig(A);
%         % Invert
%         ev = flipdim(ev,2);
%         ew = flipdim(diag(ew),1);
%         sig = ew >= podsettings.value;
%         % Reduce
%         ev = ev(:,sig);
%         ew = ew(sig);
%     % Reduce dimensions to 'value' percent
%     case 'rel'
%         % Model dimension is first dimension of snapshots
%         reddim = round(size(sn,1)*podsettings.value);
%         if reddim > size(A,2)
%             warning('KerMor:POD:mode_rel','Reduced space cant be bigger (%d) than sample size (%d). Setting to sample size.',reddim,size(A,2));
%             reddim = size(A,2);
%         end
%         % generate descending list of eigenvectors/values:
%         [ev, ew] = eigs(A,reddim);
%         ew = diag(ew);
%     % Explicit target dimension
%     case 'abs'
%         if podsettings.value > size(A,2)
%             warning('KerMor:POD:mode_abs','''model.reduction.POD.value'' cant be bigger (%d) than sample size (%d). Setting to sample size.',podsettings.value,size(A,2));
%             podsettings.value = size(A,2);
%         end
%         [ev, ew] = eigs(A,podsettings.value);
%         ew = diag(ew);
%     otherwise
%         error(['unknown POD reduction mode: ''' podsettings.mode '''']);
% end

% Compute projection matrix
model_data.V = sn * ev * diag(ew.^(-0.5));

end

