% @TODO message-system über alle berechnungen hinaus (ungewöhliche dinge
% berichten, exit flags etc)
% @TODO laufzeittests für reduzierte modelle
% @TODO interface für ModelData/Snapshots -> entweder arbeiten auf der
% Festplatte oder in 4D-Array .. (für große simulationen)
% @TODO mehr tests / anwendungen für mehrere inputs aber keine parameter!
% @TODO Verbose-Level benutzen / anpassen
% @TODO test für rotationssensitive kerne!
% @TODO: snapshotgenerierung -> mit fehlerschätzer auswählen! (3-4
% zufällig, dann approx, fehler -> neuen snapshot beim größten fehler etc)
% @TODO: moving least squares (mit gewichtsfkt) für general.regression ..
% -> book scattered data approx
% @TODO: fft-approximation (?)
% @TODO: Kern mit `\Phi(x,y) = (1-||x-y||_2)_+` oder so
% @TODO: p-partitioning
% @TODO: adaptive svr (impl. \nu-SVR, dann snapshots adden bis tol
% erreicht)
% @TODO: zusammenlegen von funktionen / erstellen eines general-modules für
% KerMor/rbmatlab?
% @TODO: try..catch langsam?
% test: zeitabhängige outputconvertierung?
% testing.MUnit auch für "nicht-packages"
% datenhaltung auf festplatte (mu,inputidx-indiziert) (?) => 
%   - berechnung kernmatrix in teilen...
%   - hashfunktion bernard / ggf eigene interface-fkt für eindeutige dirnames
% parfor für sampling / comp-wise approximation? (snaphshot-generation/approx)
 

% preferences
setpref('Internet','SMTP_Server','localhost');

% get current directory;
disp('Starting up KerMor in directory:');
p = fileparts( which('startup_kermor'));
disp(p);

% add further paths to MATLABPATH
addpath(p);

% For PCAfixspace
addpath('/afs/.mathe/project/agh/home/dwirtz/rbmatlab/general/vecmat');

cd(p);
clear('p');