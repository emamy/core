classdef KerMor < handle
    % Global configuration class for all KerMor run-time settings.
    %
    %
    % To-Do's for KerMor:
    %
    % @todo message-system über alle berechnungen hinaus (ungewöhliche dinge
    % berichten, exit flags etc)
    %
    % @todo laufzeittests für reduzierte modelle
    %
    % @todo interface für ModelData/Snapshots -> entweder arbeiten auf der
    % Festplatte oder in 4D-Array .. (für große simulationen) -> KerMor hat
    % string für globales datenverzeichnis!
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
    % @todo: datenhaltung auf festplatte (mu,inputidx-indiziert) (?) =>
    %   - berechnung kernmatrix in teilen...
    %   - hashfunktion bernard / ggf eigene interface-fkt für eindeutige dirnames
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
    % @todo speichermanagement: große matrizen / etc virtuell auf festplatte
    % laden/speichern
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
    % @todo: umstellen von simulate(mu,inputidx) auf simulate +
    % setMu,setInputidx -> faster evaluation
    
    % DONE Allgemeineres Skalarprodukt def. über `<x,y>_G = x^tGy`, default Id
    % DONE Allgemeinere Projektion mit `V,W` und nicht mit `V,V^t`
    % DONE fehler in ODE mit reinformulieren!
    % DONE getConfig-methode: string-ausgabe aller einstellungen (sofern
    % textuell sinnvoll möglich!) eines Modells
    
    properties
        % The default directory to use for simulation data storage
        DataStoreDirectory = '/datastore';
        
        % The Verbose Mode for KerMor.
        % The higher Verbose is set, the more output is produced.
        Verbose = 1;
    end
    
    methods(Static)
        function theinstance = Instance
            persistent instance;
            if isempty(instance)
                instance = KerMor;
            end
            theinstance = instance;
        end
        
        function setVerbose(value)
            k = KerMor.Instance;
            k.Verbose = value;
        end
    end
    
    methods(Access=private)
        function this = KerMor
        end
    end
    
end