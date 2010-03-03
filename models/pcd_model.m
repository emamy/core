function model = pcd_model
%BASE_MODEL The basic model every KerMor model starts from
%   Detailed explanation goes here

%% Overall model struct
model = base_model;
model.info.name = 'Programmed Cell Death';
% Determine verbose mode
model.info.verbose = 1; % Levels 0=none, 1=some, 2=extensive

%% The underlying dynamical system
model.system = pcd_dynsystem;

%% Boundary settings
model.T = 1; % Max. Simulation time
model.timestep = .1; % Timesteps to use
model.times = 0:model.timestep:model.T;

%% Sampling settings
% README: For any given total number "model.sampling.samples" the
% samples created are weighted with the corresponding
% "Desired"-fields of the param structs.
% Sampling mode: Allowed: 'grid' and 'rand'
model.sampling.mode = 'grid';
% For 'rand': Number of samples
model.sampling.samples = 10;

%% ODE Solver settings
model.odesolver = @ode45; %@ode23s; %

%% Reduction settings
% The mode determines the algorithm used to compute the projection space.
% Possible values: 'POD',
model.reduction.mode = 'POD';
%% POD settings
% README: Two possible values: 'frac' and 'abs'
% 'sign': All eigenvectors with evalues larger than 'value' percent of the
%         largest eigenvalue are used. Use 0 for all.
% 'eps':  All eigenvectors with evalues larger than 'value' are used.
% 'rel' : Reduction to 'value' percent of the original space dimension.
% 'abs' : Explicitly specified reduced space dimension. Uses the first
%         'value' eigenvectors as space.
model.reduction.POD.mode = 'eps';
model.reduction.POD.value = .1;

%% Kernel settings
% Main kernel for approximation
model.kernels.combined = @triple_product_kernel;
% Specific kernels for each variable/component
model.kernels.time = get_rbf_kernelfun(1);
model.kernels.space = get_rbf_kernelfun(1);
model.kernels.param = get_rbf_kernelfun(1);

%% Approximation settings for nonlineariy f
model.approx.mode = 'scalar_svr';
model.approx.scalar_svr.eps = .3;
model.approx.scalar_svr.C = 1;

%% Model main functions
model.plot = @plotSnapshots;
%model.gen_model_data = @generate_training_samples;


    function plotSnapshots(model, model_data)
        sns = model_data.snapshots;
        figure(1);
        d1 = model.system.d1;
        d2 = model.system.d2;
        m = d1*d2;
        for pidx = 1:size(sns,3)
            sn = sns(:,:,pidx);
            for tidx = 1:size(sn,2)
                s = sn(:,tidx);
                
                subplot(2,2,1);
                surf(reshape(s(1:m),d1,d2));
                title(sprintf('caspase-8, t=%f',model.times(tidx)));
                
                subplot(2,2,2);
                surf(reshape(s(m+1:2*m),d1,d2));
                title(sprintf('caspase-3, t=%f',model.times(tidx)));
                
                subplot(2,2,3);
                surf(reshape(s(2*m+1:3*m),d1,d2));
                title(sprintf('procaspase-8, t=%f',model.times(tidx)));
                
                subplot(2,2,4);
                surf(reshape(s(3*m+1:end),d1,d2));
                title(sprintf('procaspase-3, t=%f',model.times(tidx)));
                
                pause(1);
            end
        end
    end

end

