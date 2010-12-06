% Preferences
setpref('Internet','SMTP_Server','localhost');

% Environment
disp('Initializing environment...')
p = fileparts( which('startup_kermor'));
fprintf('Starting up KerMor in directory %s\n',p);

setenv('KERMORTEMP','/datastore');
setenv('KERMORHOME',p);

% For PCAfixspace
addpath('/afs/.mathe/project/agh/home/dwirtz/Software/rbmatlab/general/vecmat');

% add further paths to MATLABPATH
addpath(p);
addpath(fullfile(p,'3rdparty'));
addpath(fullfile(p,'3rdparty','ipopt'));
addpath(fullfile(p,'3rdparty','qpOASES'));
addpath(fullfile(p,'3rdparty','mosek'));
addpath(fullfile(p,'3rdparty','pardiso'));
addpath(fullfile(p,'demos'));
addpath(fullfile(p,'visual'));
disp('Entering startup path..');
cd(p);
clear('p');

disp('Starting up parallel computing..');
matlabpool open;
