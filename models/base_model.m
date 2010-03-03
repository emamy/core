function model = base_model
%BASE_MODEL The basic model every KerMor model starts from
%   Detailed explanation goes here

%% Overall model struct
model = struct;
model.info.name = 'KerMor base model';
% Determine verbose mode
model.info.verbose = 0; % Levels 0=none, 1=some, 2=extensive

%% The underlying dynamical system
model.system = base_dynsystem;


%% Boundary settings
model.T = .5; % Max. Simulation time
model.timestep = .1; % Timesteps to use
model.times = 0:model.timestep:model.T;

%% Sampling settings
% README: For any given total number "model.sampling.samples" the
% samples created are weighted with the corresponding
% "Desired"-fields of the param structs.
% Sampling mode: Allowed: 'grid' and 'rand'
model.sampling.mode = 'rand';
% For 'rand': Number of samples
model.sampling.samples = 10;

%% ODE Solver settings
model.odesolver = @ode45;

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
model.reduction.POD.mode = 'rel';
model.reduction.POD.value = .3;

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
%model.gen_model_data = @generate_training_samples;


end

