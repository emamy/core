classdef KerMor < handle
    % Global configuration class for all KerMor run-time settings.
    %
    % Software documentation can be found at
    % http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    %
    % Any KerMor developers should check out the pages @ref development for coding guidelines and
    % conventions.
    %
    % @author Daniel Wirtz @date 2011-03-04
    %
    % @change{0,5,dw,2011-07-07} New property DefaultFigurePosition.
    %
    % @change{0,5,dw,2011-06-20} Started KerMor version 0.5
    %
    % @change{0,4,dw,2011-05-03} 
    % - Committed KerMor version 0.3 to GIT.
    % - Different documentation creation subroutines for unix and windows now, windows not supported
    % yet
    % - Version Number now gets passed to the documentation creation script
    %
    % @new{0,3,dw,2011-04-26} Added a new property DocumentationLocation.
    %
    % @change{0,3,sa,2011-04-14} Implemented UseMatlabParallelComputing functionality
    %
    % @change{0,3,dw,2011-04-12} 
    % - Included the setup script for Devel in the installation procedure.
    % - Added new properties MainVersion and SubVersion to this class for global versioning.
    %
    % @new{0,3,dw,2011-04-05} Added a new property KerMor.DesktopLayout. As with KDE4 matlab seems
    % to start with a pretty random layout this property enables the developer to specify a
    % (previously saved) desktop layout to be associated with KerMor.
    %
    % @change{0,3,dw,2011-03-22} Modified the GIT repository so that branches now reflect the KerMor
    % versions.
    %
    % @change{0,2,dw,2011-03-21} Nicer color output for the KerMor.createDocs command and the log
    % file is printed directly in MatLab.
    %
    % @new{0,2,dw,2011-03-17} Added the fields @ref KerMor.Hasrbmatlab and KerMor.rbmatlabDirectory.
    % This allows to register a copy of rbmatlab located at the 'rbmatlabDirectory' with KerMor.
    % Without having setup this directory any models.rbmatlab classes will not work correctly.
    %
    % @change{0,2,dw,2011-03-09} Included a developer key parameter into the '@@new' and
    % '@@change{0' tags to create a link to the author of the new features or changes.
    %
    % @new{0,2,dw,2011-03-09}
    % - Added installation routines for unix systems. Now one can download the sources from a git
    % repository and simply call KerMor.install to prepare the environment & compile any included
    % mex files.
    % - Added a KerMor.createDocs static method to create the documentation from within the matlab
    % environment.
    %
    % @new{0,2,dw,2011-03-08} Created a new Doxygen keyword '@@new' for new
    % feature versioning lists
    %
    % @change{0,2,dw,2011-03-04} Moved startup and shutdown functions into
    % this class
    %
    % @change{0,1,dw} Initial version. @new{0,1,dw} Initial version.
    %
    % To-Do's for KerMor:
    % @todo
    % - message-system ?ber alle berechnungen hinaus (ungew?hliche dinge
    % berichten, exit flags etc). hier eine zentrale logging-funktion die
    % je nach verbose sachen direkt plottet oder (immer!) in ein log-file
    % schreibt.. das ganze mit verbose-leveln kombinieren!
    % - laufzeittests f?r reduzierte modelle
    % - interface f?r ModelData/Snapshots -> entweder arbeiten auf der
    % Festplatte oder in 4D-Array .. (f?r gro?e simulationen) -> KerMor hat
    % string f?r globales datenverzeichnis!
    % generell datenhaltung auf festplatte (mu,inputidx-indiziert) (?) =>
    %   - berechnung kernmatrix in teilen...
    %   - hashfunktion bernard / ggf eigene interface-fkt f?r eindeutige
    %   dirnames
    %   -speichermanagement: gro?e matrizen / etc virtuell auf festplatte
    %   laden/speichern
    %
    % @todo mehr tests / anwendungen f?r mehrere inputs aber keine parameter!
    %
    % @todo Verbose-Level benutzen / anpassen
    %
    % @todo test f?r rotationssensitive kerne!
    %
    % @todo: snapshotgenerierung -> mit fehlersch?tzer ausw?hlen! (3-4
    % zuf?llig, dann approx, fehler -> neuen snapshot beim gr??ten fehler etc)
    %
    % @todo: moving least squares (mit gewichtsfkt) f?r general.regression ..
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
    % @todo: zusammenlegen von funktionen / erstellen eines general-modules f?r
    % KerMor/rbmatlab?
    %
    % @todo: try..catch langsam?
    %
    % @todo zeitabh?ngige outputconvertierung?
    % testing.MUnit auch f?r "nicht-packages"
    %
    %
    % @todo: parfor f?r sampling / comp-wise approximation? (snaphshot-generation/approx)
    %
    % @todo benchmarks von
    % http://portal.uni-freiburg.de/imteksimulation/downloads/benchmark
    % einlesbar machen / einbauen!
    %
    % @todo Beispiele von ODE's aus Matlab-Docs? (verficiation)
    %
    % @todo Fehlersch?tzer auf Output beschr?nken/erweitern!
    %
    % @todo Mehr ODE-Solver (implizit) einbauen, ggf. eigenen RK23 oder so.
    %
    % @todo LaGrange-koeffizientenfunktionen bei kerninterpolation berechnen!
    % ist insgesamt billiger falls `N<<n`
    %
    % @todo: test f?r newton-iteration!
    %
    % @todo Implementierung Balanced Truncation (mit base class) f?r
    % LinearCoreFuns, dann implementierung balanced truncation f?r empirical
    % gramians nach paper Lall et al. -> neue subspace reduction method f?r
    % nonlin-systems mit inputs! (geht ggf. auch f?r systeme ohne inputs? probieren!)
    %
    % @todo vielleicht so etwas wie "isValid" f?r jede modellkomponente, das
    % vor start von teuren berechnungen pr?ft ob alles so durchgeht und keine
    % inkompatibilit?ten auftreten (z.B. krylov - LinearCoreFun)
    %
    % @todo check ob es eine m?glichkeit gibt zu pr?fen ob alle unterklassen
    % von projizierbaren klassen die project-methode der oberklasse aufrufen!?
    % k?nnte sonst zu komischen fehlern f?hren..
    %
    % @todo check warum der error estimator nach einem save vom reduzierten
    % modell nicht gespeichert wird.
    %
    % @todo 16.09.2010: f?r skalarprodukt-kerne eigenes interface
    % implementieren und dann ggf. f?r W=V korrekt projezieren + TEST
    % schreiben!
    %
    % @todo cacher aus RBMatlab portieren/?bertragen!
    %
    % @todo t-partitioning f?r KerMor? ideen mit markus austauschen!
    %
    % @todo MUnit erweitern um benchmark-mode, so das anstelle von "test_"
    % prefix-fkt alle mit "bench_" ausgef?hrt werden; (r?ckgabe ist in dem fall
    % ggf ein struct mit algorithmus und zeiten)
    %
    % @todo eigene POD-Basen f?r verschiedene Teile des systems denen andere
    % physik zugrunde liegt (i.e. f(x) => [f_1(x); f_2(x)]), mit einzelner
    % auswertung? dazu m?sste man indexmatrizen einrichten die die
    % verschiedenen teile von f bezeichnen... (Motivation: "Application of POD
    % and DEIM for MOR of Nonl. Miscible Viscous Flow, Chaturantabut/Sorensen)
    %
    % @todo fehlersch?tzer gegen die volle, nicht projizierte
    % kernelapproximation einrichten? damit kann man den aktuell besch?tzten
    % fehler besser bekommen..
    %
    % @todo sekantenabsch?tzung per kernregression vorab f?r 1D berechnen? dann
    % entf?llt das newton-problem l?sen. interpolation z.B. geht auch f?r
    % fehlerabsch?tzung um die rigorosit?t zu erhalten.
    %
    % @todo timedirty ?berarbeiten / rausnehmen etc, sollte auch einzelaufrufe
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
    
    properties(Constant)
        % The current KerMor main version number
        %
        % Change only AFTER committing the final last version's state.
        % Used in Devel to fill the new class templates etc.
        %
        % See also: SubVersion
        MainVersion = '0';
        
        % The current KerMor sub version number
        %
        % Change only AFTER committing the final last version's state.
        % Used in Devel to fill the new class templates etc.
        %
        % See also: MainVersion
        SubVersion = '5';
    end
    
    properties
        % The directory to use for simulation data storage
        %
        % In this folder large simulation and model data will be stored and
        % read by KerMor. Can be set anytime during runtime. If no value is
        % given the KerMor.start script will ask for it.
        %
        % @default ./data
        DataStoreDirectory = '';
        
        % The directory to use for temporary simulation data
        % @default ./temp
        TempDirectory = '';
        
        % The preferred desktop layout to work with.
        %
        % If you work with different desktop layouts or the KDE JUST DOES
        % NOT GET IT you can save your custom desktop layout and set this
        % property to its name. Upon start, KerMor will restore the layout
        % for you automatically. Set to '' to disable.
        % @default empty
        DesktopLayout = '';
        
        % The source directory for a copy of rbmatlab
        %
        % @default []
        rbmatlabDirectory = '';
        
        % Verbose output level
        %
        % @default 1
        Verbose = 1;
        
        % Flag whether to enable use of the Matlab Parallel Computing
        % Toolbox.
        %
        % @default false
        UseMatlabParallelComputing = false;
        
        % The default figure position to use.
        %
        % If none is set, KerMor does not modify the root workspace property
        % 'DefaultFigurePosition' upon startup.
        DefaultFigurePosition = [];
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
        
        % Flag if rbmatlab wrapping functionalities are enabled
        %
        % @default false
        Hasrbmatlab = false;
    end
    
    properties(Dependent)
        % Returns where the documentation is located.
        %
        % This is either a web-site or the local documentation.
        DocumentationLocation;
    end
    
    methods
        function set.UseMatlabParallelComputing(this, value)
            if ~islogical(value)
                error('Value must be logical');
            end
            haspc = ~isempty(which('matlabpool'));
            if haspc
                s = matlabpool('size');
            end
            if value            
                if haspc
                   this.UseMatlabParallelComputing = value;
                   setpref('KERMOR','USEMATLABPARALLELCOMPUTING',value);
                   if s == 0
                       matlabpool open;
                   end                
                else
                    error('No parallel computing toolbox available.');
                end
            else
                this.UseMatlabParallelComputing = value;
                setpref('KERMOR','USEMATLABPARALLELCOMPUTING',value);
                if s > 0
                   matlabpool close;
                end
            end
        end
          
        function set.DataStoreDirectory(this, value)
            if ~isempty(value) && ~isdir(value)
                fprintf('Creating directory %s\n',value);
                mkdir(value);
            end
            setpref('KERMOR','DATASTORE',value);
            this.DataStoreDirectory = value;
            fprintf('Simulation and model data: %s\n',value);
        end
        
        function set.TempDirectory(this, value)
            if ~isempty(value) && ~isdir(value)
                fprintf('Creating directory %s\n',value);
                mkdir(value);
            end
            setpref('KERMOR','TMPDIR',value);
            this.TempDirectory = value;
            fprintf('Temporary files: %s\n',value);
        end
        
        function set.DesktopLayout(this, value)
            setpref('KERMOR','DESKLAYOUT',value);
            this.DesktopLayout = value;
        end
        
        function set.rbmatlabDirectory(this, value)
            % Sets the rbmatlab source directory
            %
            % Parameters:
            % ds: The directory of an rbmatlab source. Use empty string or
            % cell to "uninstall" rbmatlab.
            %
            % Throws an exception if the path is invalid or does not
            % contain the rbmatlab startup script.
            if ~isempty(value)
                if ~isdir(value)
                    error('Invalid directory: %s',value);
                end
                if ~exist(fullfile(value,'startup_rbmatlab.m'),'file')
                    error('Invalid rbmatlab directory (no startup script found): %s',value);
                end
            end
            setpref('KERMOR','RBMATLABDIR',value);
            this.rbmatlabDirectory = value;
            fprintf('rbmatlab root directory: %s\n',value);
        end
        
        function value = get.UseMatlabParallelComputing(this)
            % recover values if clear classes has been issued or Matlab
            % Parallel processing toolbox deleted            
            value = getpref('KERMOR','USEMATLABPARALLELCOMPUTING','');
            t = which('matlabpool');
            if ~isempty(value)
                if ~isempty(t)
                    this.UseMatlabParallelComputing = value;
                elseif value
                    warning('KERMOR:ParComp','No parallel computing toolbox found but preference for UseMatlabParallelComputing contained true. Setting to false.');
                    this.UseMatlabParallelComputing = false;
                end
            else
                error('Invalid configuration detected. Have you run the installation routine?');
            end
        end
        
        function h = get.HomeDirectory(this)
            if isempty(this.HomeDirectory)
                this.HomeDirectory = fileparts(which('KerMor'));
            end
            h = this.HomeDirectory;
        end
        
        function d = get.DocumentationLocation(this)%#ok
            d = getenv('KERMOR_DOCS');
            if isempty(d) || ~exist(fullfile(d,'index.html'),'file')
                d = 'http://www.agh.ians.uni-stuttgart.de/documentation/kermor';
            end
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
        
        function d = get.DesktopLayout(this)
            % recover values if clear classes has been issued
            if isempty(this.DesktopLayout)
                d = getpref('KERMOR','DESKLAYOUT','');
                this.DesktopLayout = d;
            else
                d = this.DesktopLayout;
            end
        end
        
        function h = get.rbmatlabDirectory(this)
            
            % recover values if clear classes has been issued
            if isempty(this.rbmatlabDirectory)
                h = getpref('KERMOR','RBMATLABDIR','');
                if ~isempty(h)
                    this.rbmatlabDirectory = h;
                end
            else
                h = this.rbmatlabDirectory;
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
        
        function flag = get.Hasrbmatlab(this)
            % Indicates if rbmatlab is available and initialized
            %
            % This function checks for a set rbmatlabDirectory and if the
            % rbmatlab script 'rbmatlabhome' is within the current path
            % (=substitute for the rbmatlab startup-script being executed)
            flag = false;
            if ~isempty(this.rbmatlabDirectory)
                flag = ~isempty(which('rbmatlabhome'));
                if ~flag
                    warning('KerMor:App',...
                        ['rbmatlab directory is set, but script'...
                        ' ''rbmatlabhome'' could not be found in current path.\n'...
                        'Unsure if rbmatlab-dependent models will work,'...
                        ' check if rbmatlab version has changed!']);
                end
            end
        end
        
        function value = get.DefaultFigurePosition(this)
            value = this.DefaultFigurePosition;
            if isempty(value)
                value = getpref('KERMOR','DefFigPos',[]);
                if ~isempty(value)
                    this.DefaultFigurePosition = value;
                end
            end
        end
        
        function set.DefaultFigurePosition(this, value)
            setpref('KERMOR','DefFigPos',value);
            this.DefaultFigurePosition = value;
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
            addpath(fullfile(p,'external'));
            addpath(fullfile(p,'external','WH10'));
            addpath(fullfile(p,'external','ICIAM2011'));
            % Figure position settings
            if ~isempty(this.DefaultFigurePosition)
                set(0,'DefaultFigurePosition',this.DefaultFigurePosition);
            end
            
            initDirectories;
            init3rdparty;            
            initParallelization;
            
            disp('Entering startup path..');
            cd(p);
            clear('p');
            
            if ~isempty(this.DesktopLayout)
                fprintf('Applying desktop layout %s..\n',this.DesktopLayout);
                desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
                desktop.restoreLayout(this.DesktopLayout);
            end
            
            
            
            
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
                
                % md5
                %addpath(fullfile(p,'3rdparty','md5'));
                
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
                
                % rbmatlab
                if ~isempty(this.rbmatlabDirectory)
                    disp('<<<<<<<<< Starting rbmatlab >>>>>>>>>>');
                    setenv('RBMATLABTEMP', this.TempDirectory);
                    setenv('RBMATLABHOME', this.rbmatlabDirectory);
                    addpath(this.rbmatlabDirectory);
                    curdir = pwd;
                    evalin('base','startup_rbmatlab;');
                    evalin('base','clear all;');
                    chdir(curdir);
                    disp('<<<<<<<<< Done Starting rbmatlab >>>>>>>>>>');
                end
                
                %addpath(fullfile(p,'3rdparty','pardiso'));
            end   
            
            function initParallelization
                % Checks if the parallel computing toolbox is available
                %
                % @todo wrap with try-catch and set flag in KerMor.App
                % class!
                %
                % @note The 'feature' command is undocumented, see
                % http://www.mathworks.com/matlabcentral/newsreader/view_thread/154551
                % for more information.
                disp('Checking for and starting parallel computing..');
                
                % Open matlabpool only if UseMatlabParallelComputing is set to
                % true
                if this.UseMatlabParallelComputing
                    if matlabpool('size') == 0
                        matlabpool open;
                    end
                end
                
                % Sets the maximum number of threads to create by OpenMP
                % binaries according to the number of cores available on
                % the machine.
                setenv('OMP_NUM_THREADS',num2str(feature('numCores')));
            end
        end        
                  
        function shutdown(this)%#ok
            % only close if parallel computing is available and matlabpool is running!
            t = which('matlabpool');
            if ~isempty(t)
                s = matlabpool('size');
                if s > 0
                    matlabpool close;
                end
            end           
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
        
        function install
            % Performs installation of KerMor on a system
            %
            % Adds variables to the users environment, so far only needed
            % for the documentation creation. Custom paths for data storage
            % are checked and set by the start script.
            %
            % If no rbmatlab directory is already presently set by a
            % previous install, the installation program asks if an
            % rbmatlab-installation should be registered with KerMor.
            %
            % See also: installUnix installWindows
            disp('<<<<<<<<<< Welcome to the KerMor install script. >>>>>>>>>>');
            
            %% Operation-system dependent actions
            if isunix
                KerMor.installUnix;
            elseif ispc
                KerMor.installWindows;
            end
            
            %% Setup KerMor development
            if isempty(getpref('KERMOR_DEVEL','author',''))
                str = sprintf('Do you want to setup variables for KerMor development?\n(Y)es/(N)o: ');
                ds = lower(input(str,'s'));
                if isequal(ds,'y')
                    Devel.setup;
                end
            end
            
            %% SETPREF for Matlab Parallel processing
            if ~isempty(which('matlabpool'));
                str = sprintf('Do you want to Use Matlab Parallel Processing?\n(Y)es/(N)o: ');
                value = lower(input(str,'s'));
                if isequal(value,'y')
                    setpref('KERMOR','USEMATLABPARALLELCOMPUTING','true');
                else
                    setpref('KERMOR','USEMATLABPARALLELCOMPUTING','false');
                end
            end
            %% Optional: rbmatlab
            a = KerMor.App;
            if isempty(a.rbmatlabDirectory)
                str = sprintf(['Do you want to register a local rbmatlab '...
                    ' version with KerMor?\n(Y)es/(N)o: ']);
                ds = lower(input(str,'s'));
                if isequal(ds,'y')
                    d = uigetdir(h,'Please select the rbmatlab source root folder.');
                    if d ~= 0
                        try
                            a.rbmatlabDirectory = d;
                        catch ME
                            disp('Setting the rbmatlab directory failed:');
                            disp(getReport(ME, 'basic'));
                            disp('You can still connect to rbmatlab later by setting the rbmatlabDirectory property.');
                        end
                    end
                end
            end
            disp('<<<<<<<<<< Setup complete. You can now start KerMor by running "KerMor.start;". >>>>>>>>>>');
        end
        
        function createDocs(uml, open)
            % Creates the Doxygen documentation
            %
            % Parameters:
            % uml: Set to true to create UML-like graphics output
            % @default false
            % open: Set to true if the documentation should be opened after
            % successful compilation
            % @default false
            if nargin < 2
                open = false;
                if nargin < 1
                    uml = false;
                end
            end
            
            %% Operation-system dependent actions
            if isunix
                KerMor.createDocsUnix(uml, open);
            elseif ispc
                KerMor.createDocsWindows(uml, open);
            end
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
    
    methods(Static,Access=private)
        
        function createDocsUnix(uml, open)
            cmd = fullfile(KerMor.App.HomeDirectory,'documentation','make.sh');
            % Add version argument
            cmd = [cmd ' ' KerMor.MainVersion '.' KerMor.SubVersion];
            % Add uml argument
            if uml
                cmd = [cmd ' uml'];
            end
            [s,r] = system(cmd);
            wpos = strfind(r,'Logged warnings:');
            if ~isempty(wpos)
                endpos = strfind(r,'Complete log file');
                cprintf([0 .5 0],r(1:wpos-1));
                cprintf([1,.4,0],strrep(r(wpos:endpos-1),'\','\\')); 
                cprintf([0 .5 0],r(endpos:end));
            else
                cprintf([0 .5 0],strrep(r,'\','\\'));
            end
            fprintf('\n');
            index = fullfile(getenv('KERMOR_DOCS'), 'index.html');
            if open
                % Try to use iceweasel per default (nasty, uhm?)
                [s,r] = system('which iceweasel');
                if ~isempty(r)
                    cmd = 'iceweasel ';
                else
                    % Otherwise: use user preferred browser
                    cmd = 'xdg-open ';
                end
                system([cmd index]);
            end
        end
        
        function createDocsWindows(uml, open)%#ok
            error('Creating documentation on Windows is currently not supported.');
        end
        
        function installUnix
            % Install script for unix systems
            %
            % Adds kermor specific variables to the user's environment by
            % inserting them into the ~/.bashrc file.
            %
            % @note If you run this install script more than once, old path
            % variables will be overwritten (as multiple entries will
            % appear in the .bashrc file)
            %
            % The custom variable names are
            % - 'KERMOR_SOURCE' The source directory
            % - 'KERMOR_DOCS' The documentation output directory
            % - 'KERMOR_DOXYBIN' Path to the doxygen binary (autodetect)
            %
            % @todo compile any mex files!
            
            h = fileparts(which('KerMor'));
            fid = fopen('~/.bashrc','a+');
            try
                
                fprintf(fid,'\n# KerMor environment variables (added by KerMor.install script on %s)\n',date);
                fprintf(fid,'export KERMOR_SOURCE="%s"\n',h);
                % Set in running environment (until restart)
                setenv('KERMOR_SOURCE',h);
                
                
                %% Documentation directory
                if isempty(getenv('KERMOR_DOCS'))
                    d = fullfile(h,'documentation','output');
                    str = sprintf(['No documentation output directory has been set yet.\n'...
                        'The default will be %s\n'...
                        'Do you want to specify a custom output directory? (Y)es/(N)o: '],d);
                    ds = lower(input(str,'s'));
                    if isequal(ds,'y')
                        d = uigetdir(h,'Please select the documentation output folder.');
                        if d == 0
                            d = fullfile(h,'documentation','output');
                            fprintf('Operation cancelled, using default directory %s...\n',d);
                        end
                    end
                    fprintf(fid,'export KERMOR_DOCS="%s"\n',d);
                    setenv('KERMOR_DOCS',d)
                end
                
                %% Doxygen binary
                if isempty(getenv('KERMOR_DOXYBIN'))
                    [s,r] = system('which doxygen');
                    db = [];
                    if ~isempty(r)
                        db = 'doxygen';
                        str = sprintf(['Doxygen installation is available (%s).\n'...
                            'Do you want to use a different doxygen binary? (Y)es/(N)o: '],strrep(r,char(10),''));
                        yn = lower(input(str,'s'));
                    end
                    if isempty(r) || isequal(yn,'y')
                        [f,p] = uigetfile('~/*.*','Select the custom doxygen binary file.');
                        if f ~= 0
                            db = fullfile(p, f);
                        end
                    end
                    if ~isempty(db)
                        fprintf(fid,'export KERMOR_DOXYBIN="%s"\n',db);
                        setenv('KERMOR_DOXYBIN',db);
                    else
                        warning('KerMor:installUnix','No doxygen binary selected. Documentation creation will not work.');
                    end
                end
                
                fclose(fid);
                
                % Reload terminal bashrc
                system('. ~/.bashrc');
            catch ME%#ok
                fclose(fid);
            end
        end
        
        function installWindows
            % installation script for Microsoft Windows based systems.
            %
            % Not yet implemented/necessary (docs make-script is linux only!)
            %
            % @todo
            % - Install script for Windows
            % - Create batch file for documentation creation on windows
            error('Installation routine not yet implemented.\nPlease refer to the KerMor documentation at http://www.agh.ians.uni-stuttgart.de/documentation/kermor for help.');
        end
    end
    
    methods(Access=private)
        function this = KerMor
            % Private constructor: This class is a Singleton.
        end
    end
    
end