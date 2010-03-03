% preferences
setpref('Internet','SMTP_Server','localhost');

% get current directory;
disp('Starting up KerMor in directory:');
p = fileparts( which('startup_kermor'));
disp(p);

% add further paths to MATLABPATH
addpath( fullfile( p ,'general') );
addpath( fullfile( p ,'kernels') );
addpath( fullfile( p ,'models') );
addpath( fullfile( p ,'svr') );
addpath( fullfile( p ,'testing') );

% For PCAfixspace
addpath('/afs/.mathe/project/agh/home/dwirtz/rbmatlab/general/vecmat');

cd(p);
clear('p');