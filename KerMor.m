classdef KerMor < handle
    % Global configuration class for all KerMor run-time settings.
    %
    % Software documentation can be found at
    % http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    %
    % @date 04.03.2011 @author Daniel Wirtz
    %
    % @change{0,2,2011-03-04} Moved startup and shutdown functions into this class
    % @change{0,2,2011-03-03} Created local settings file and management thereof for
    % easy local installation and path adjustment.
    % @change{0,1} Initial version.
    %
    % To-Do's for KerMor:
    % @todo message-system über alle berechnungen hinaus (ungewöhliche dinge
    % berichten, exit flags etc)
    %
    % @todo laufzeittests für reduzierte modelle
    %
    % @todo interface für ModelData/Snapshots -> entweder arbeiten auf der
    % Festplatte oder in 4D-Array .. (für große simulationen) -> KerMor hat
    % string für globales datenverzeichnis!
    % generell datenhaltung auf festplatte (mu,inputidx-indiziert) (?) =>
    %   - berechnung kernmatrix in teilen...
    %   - hashfunktion bernard / ggf eigene interface-fkt für eindeutige
    %   dirnames
    %   -speichermanagement: große matrizen / etc virtuell auf festplatte
    %   laden/speichern
    %
    % @todo mehr tests / anwendungen für mehrere inputs aber keine parameter!
    %
    % @todo Verbose-Level benutzen / anpassen
    %
    % @todo test für rotationssensitive kerne!
    %
    % @todo: snapshotgenerierung -> mit fehlerschätzer auswählen! (3-4
    % zufällig, dann approx, fehler -> neuen snapshot beim größten fehler etc)
    %
    % @todo: moving least squares (mit gewichtsfkt) für general.regression ..
    % -> book scattered data approx
    %
    % @todo: fft-approximation (?)
    %
    % @todo: Kern mit `\Phi(x,y) = (1-||x-y||_2)_+` oder so
    %
    % @todo: p-partitioning
    %
    % @todo: adaptive svr (impl. `\nu`-SVR, dann snapshots adden bis tol
    % erreicht)
    %
    % @todo: zusammenlegen von funktionen / erstellen eines general-modules für
    % KerMor/rbmatlab?
    %
    % @todo: try..catch langsam?
    %
    % @todo zeitabhängige outputconvertierung?
    % testing.MUnit auch für "nicht-packages"
    %
    %
    % @todo: parfor für sampling / comp-wise approximation? (snaphshot-generation/approx)
    %
    % @todo benchmarks von
    % http://portal.uni-freiburg.de/imteksimulation/downloads/benchmark
    % einlesbar machen / einbauen!
    %
    % @todo Beispiele von ODE's aus Matlab-Docs? (verficiation)
    %
    % @todo Fehlerschätzer auf Output beschränken/erweitern!
    %
    % @todo Mehr ODE-Solver (implizit) einbauen, ggf. eigenen RK23 oder so.
    %
    % @todo LaGrange-koeffizientenfunktionen bei kerninterpolation berechnen!
    % ist insgesamt billiger falls `N<<n`
    %
    % @todo: test für newton-iteration!
    %
    % @todo Implementierung Balanced Truncation (mit base class) für
    % LinearCoreFuns, dann implementierung balanced truncation für empirical
    % gramians nach paper Lall et al. -> neue subspace reduction method für
    % nonlin-systems mit inputs! (geht ggf. auch für systeme ohne inputs? probieren!)
    %
    % @todo vielleicht so etwas wie "isValid" für jede modellkomponente, das
    % vor start von teuren berechnungen prüft ob alles so durchgeht und keine
    % inkompatibilitäten auftreten (z.B. krylov - LinearCoreFun)
    %
    % @todo check ob es eine möglichkeit gibt zu prüfen ob alle unterklassen
    % von projizierbaren klassen die project-methode der oberklasse aufrufen!?
    % könnte sonst zu komischen fehlern führen..
    %
    % @todo check warum der error estimator nach einem save vom reduzierten
    % modell nicht gespeichert wird.
    %
    % @todo 16.09.2010: für skalarprodukt-kerne eigenes interface
    % implementieren und dann ggf. für W=V korrekt projezieren + TEST
    % schreiben!
    %
    % @todo cacher aus RBMatlab portieren/übertragen!
    %
    % @todo t-partitioning für KerMor? ideen mit markus austauschen!
    %
    % @todo check ob die Norm von kernexpansionen mit offset-term b ähnlich
    % bleibt!
    %
    % @todo MUnit erweitern um benchmark-mode, so das anstelle von "test_"
    % prefix-fkt alle mit "bench_" ausgeführt werden; (rückgabe ist in dem fall
    % ggf ein struct mit algorithmus und zeiten)
    %
    % @todo eigene POD-Basen für verschiedene Teile des systems denen andere
    % physik zugrunde liegt (i.e. f(x) => [f_1(x); f_2(x)]), mit einzelner
    % auswertung? dazu müsste man indexmatrizen einrichten die die
    % verschiedenen teile von f bezeichnen... (Motivation: "Application of POD
    % and DEIM for MOR of Nonl. Miscible Viscous Flow, Chaturantabut/Sorensen)
    %
    % @todo fehlerschätzer gegen die volle, nicht projizierte
    % kernelapproximation einrichten? damit kann man den aktuell beschätzten
    % fehler besser bekommen..
    %
    % @todo sekantenabschätzung per kernregression vorab für 1D berechnen? dann
    % entfällt das newton-problem lösen. interpolation z.B. geht auch für
    % fehlerabschätzung um die rigorosität zu erhalten.
    %
    % @todo timedirty überarbeiten / rausnehmen etc, sollte auch einzelaufrufe
    % zu offX checken.
    %
    % @todo umstellen von simulate(mu,inputidx) auf simulate +
    % setMu,setInputidx -> faster evaluation
    % -gegenargument: schlechte parallelisierbarkeit bei zentralem mu/inidx
    %
    % @todo PCAFixspace wieder einbauen, um greedy-basisgen zu erlauben (->
    % generell: greedy-unterraumalgorithmus einbauen)
    %
    % @todo add tests for models with C output and custom G for 1-n
    % dimensions!
    %
    % @todo hierarchical subspace selection? -> bad that approximation is
    % trained on largest subspace, so would have to have more
    % approximations. or any possibility to project approximation into
    % sub-subspace??
    %
    % @todo !!!!!!!!!!!!!!!!!! include kernel expansion offset 'b' into
    % error estimators! (yet only valid for b=0)
    
    properties
        % The directory to use for simulation data storage
        %
        % In this folder large simulation and model data will be stored and
        % read by KerMor.
        %
        % @default ./data
        DataStoreDirectory = '';
        
        % The directory to use for temporary simulation data
        % @default ./temp
        TempDirectory = '';
        
        % Verbose output level
        %
        % @default 1
        Verbose = 1;
    end
    
    properties(SetAccess=private)
        % The KerMor home directory
        HomeDirectory;
    end
    
    properties(SetAccess=private,Dependent)
        % Flag if 3rd party IPOPT is available
        %
        % @default false
        HasIPOPT = false;
        
        % Flag if 3rd party qpOASES is available
        %
        % @default false
        HasqpOASES = false;
        
        % Flag if 3rd party qpMosek is available
        %
        % @default false
        HasqpMosek = false;
    end
    
    methods
        function set.DataStoreDirectory(this, ds)
            if ~isdir(ds)
                fprintf('Creating directory %s\n',ds);
                mkdir(ds);
            end
            setpref('KERMOR','DATASTORE',ds);
            this.DataStoreDirectory = ds;
            fprintf('Simulation and model data: %s\n',ds);
        end
        
        function set.TempDirectory(this, tmp)
            if ~isdir(tmp)
                fprintf('Creating directory %s\n',tmp);
                mkdir(tmp);
            end
            setpref('KERMOR','TMPDIR',tmp);
            this.TempDirectory = tmp;
            fprintf('Temporary files: %s\n',tmp);
        end
        
        function h = get.HomeDirectory(this)
            if isempty(this.HomeDirectory)
                this.HomeDirectory = fileparts(which('KerMor'));
            end
            h = this.HomeDirectory;
        end
        
        function h = get.DataStoreDirectory(this)
            
            % recover values if clear classes has been issued
            if isempty(this.DataStoreDirectory)
                h = getpref('KERMOR','DATASTORE','');
                if ~isempty(h)
                    this.DataStoreDirectory = h;
                end
            else
                h = this.DataStoreDirectory;
            end
        end
        
        function h = get.TempDirectory(this)
            
            % recover values if clear classes has been issued
            if isempty(this.TempDirectory)
                h = getpref('KERMOR','TMPDIR','');
                if ~isempty(h)
                    this.TempDirectory = h;
                end
            else
                h = this.TempDirectory;
            end
        end
        
        function flag = get.HasIPOPT(this)%#ok
            flag = ~isempty(which('ipopt'));
        end
        
        function flag = get.HasqpOASES(this)%#ok
            flag = ~isempty(which('qpOASES'));
        end
        
        function flag = get.HasqpMosek(this)%#ok
            flag = ~isempty(which('mosekopt'));
        end
    end
    
    methods(Access=private)
        function initialize(this)
            % Internal main startup script.
            
            disp('<<<<<<<<< Welcome to KerMor! >>>>>>>>>>');
            
            disp('Initializing environment...')
            % Preferences & Environment
            setpref('Internet','SMTP_Server','localhost');
            % Setup home directory & paths
            p = this.HomeDirectory;
            addpath(p);
            addpath(fullfile(p,'demos'));
            addpath(fullfile(p,'visual'));
            
            initDirectories;
            init3rdparty;
            initParallelToolbox;
            
            disp('Entering startup path..');
            cd(p);
            clear('p');
            
            
            disp('<<<<<<<<< Ready to go. >>>>>>>>>>');
            
            function initDirectories
                % Setup the data storage directory
                ds = this.DataStoreDirectory;
                if isempty(ds)
                    ds = input(['Please specify the KerMor data file directory.\n'...
                        'Leaving it empty creates a new "data" folder within the KerMor home directory.\n'...
                        'Absolute path: '],'s');
                    if isempty(ds)
                        ds = fullfile(this.HomeDirectory,'data');
                    end
                    this.DataStoreDirectory = ds;
                end
                
                % Setup the data storage directory
                tmp = this.TempDirectory;
                if isempty(tmp)
                    tmp = input(['Please specify the KerMor temporary data file directory.\n'...
                        'Leaving it empty creates a new "temp" folder within the KerMor home directory.\n'...
                        'Absolute path: '],'s');
                    if isempty(tmp)
                        tmp = fullfile(this.HomeDirectory,'temp');
                    end
                    this.TempDirectory = tmp;
                end
            end
            
            function init3rdparty
                % Checks for 3rd party software availability
                %
                % @todo include checks for pardiso once pardiso solver is
                % implemented/wrapped
                disp('Checking for 3rd party software...')
                addpath(fullfile(p,'3rdparty'));
                
                % ipopt
                addpath(fullfile(p,'3rdparty','ipopt'));
                if ~this.HasIPOPT
                    warning('KerMor:init','No IPOPT available!');
                    %rmpath(fullfile(p,'3rdparty','ipopt'));
                end
                
                % qpoases
                addpath(fullfile(p,'3rdparty','qpOASES'));
                if ~this.HasqpOASES
                    warning('KerMor:init','No qpOASES available!');
                    %rmpath(fullfile(p,'3rdparty','qpOASES'));
                end
                
                % mosek
                addpath(fullfile(p,'3rdparty','mosek'));
                if ~this.HasqpMosek
                    warning('KerMor:init','No mosek available!');
                    %rmpath(fullfile(p,'3rdparty','mosek'));
                end
                
                %addpath(fullfile(p,'3rdparty','pardiso'));
            end
            
            function initParallelToolbox
                % Checks if the parallel computing toolbox is available
                disp('Checking for and starting parallel computing..');
                % @todo wrap with try-catch and set flag in KerMor.App class!
                matlabpool open;
            end
        end
        
        function shutdown(this)
            % @todo only close if parallel computing is available!
            matlabpool close;
        end
    end
    
    methods(Static)
        function theinstance = App
            % The singleton KerMor instance
            %
            % Access to the main programs instance via KerMor.App!
            persistent instance;
            if isempty(instance)
                instance = KerMor;
            end
            theinstance = instance;
        end
        
        function application = start
            % Starts the KerMor application
            %
            % This static method initializes the environment and performs
            % initial availability checks for e.g. 3rd party programs and
            % matlab toolboxes.
            %
            % Additionally, some path variables are tried to be read from
            % environment variables if set. If no settings are made yet the
            % user is prompted to select them!
            %
            % Return values:
            % application: The KerMor instance.
            application = KerMor.App;
            application.initialize;
%             if ~application.Started
%                 application.initialize;
%             else
%                 error('The application is already started.');
%             end
        end
        
        function stop
            % Ends the KerMor application
            %
            % Stores some global property values in environment variables
            % (paths) or text files (misc settings)
            % Closes the matlab-pool for parallel computing is closed
            % if in use.
            KerMor.App.shutdown;
        end
    end
    
    methods(Access=private)
        function this = KerMor
            % Private constructor: This class is a Singleton.
        end
    end
    
end