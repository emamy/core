function model_data = gen_approximation( model, model_data )
%GEN_APPROXIMATION Summary of this function goes here
%   Detailed explanation goes here

xi = model_data.snapshots(:,:);
fxi = model_data.f_values(:,:);

vect = compile_triplevect(model, model_data);

% Compute kernel matrix (only to be computed once for all regressions)
kernel_matrix = model.kernels.combined(vect,vect,size(xi,1),...
                    model.kernels.time,...
                    model.kernels.space,...
                    model.kernels.param);

% Approximate each f-dimension
fdims = size(fxi,1);

% TODO: rough estimation of needed
approx_data = struct('ai',{},'b',{},'svidx',{});
svidxsum = [];
try
    wh = waitbar(0,'Initializing SVR ...');
    for idx = 1:fdims
        waitbar(idx/fdims,wh,sprintf('Performing SVR for dimension %d/%d ... %2.0f %%',idx,fdims,(idx/fdims)*100));
        
        [ai,b,svidx] = scalar_svr(fxi(idx,:), kernel_matrix,...
                                  model.approx.scalar_svr.eps,...
                                  model.approx.scalar_svr.C);
        
        approx_data(idx).ai = ai;
        approx_data(idx).b = b;
        approx_data(idx).svidx = svidx;
        % collect all used support vector indices
        svidxsum = union(svidxsum,svidx);
    end
catch ME
    close(wh);
    rethrow(ME);
end

waitbar(1,wh,'Compiling approximation model data...');

vect = vect(svidxsum);
% Create transition matrix for index updates
trans(svidxsum) = 1:length(svidxsum);
% Update the support vector indices to the new index in the reduced support
% vector set. Unfortunately can first be done after finishing scalar
% approximation (hence second for-loop)
for idx = 1:fdims
    approx_data(idx).svidx = trans(approx_data(idx).svidx);
end
model_data.suppvect = vect;
model_data.approx_data = approx_data;

if ishandle(wh)
    close(wh); 
end
                
    
end

