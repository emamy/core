classdef MLode15i < solvers.ode.BaseImplSolver
% MLode15i: Wrapper for MatLab's ode15i builtin implicit solver
%
% @author Daniel Wirtz @date 2011-04-14
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @new{0,3,dw,2011-04-14} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetObservable)
        % Relative error tolerance for solver.
        %
        % Set to a lower value than the default 1e-3 of matlab ode solvers.
        %
        % See the RelTol description in the odeset documentation for the effects of this property.
        %
        % @propclass{critical} The correct relative tolerance guarantees the right number of
        % significant figures.
        %
        % @default 1e-4;
        %
        % See also: AbsTol ode15i odeset
        RelTol = 1e-4;
        
        % Absolute error tolerance for solver.
        %
        % Set to a lower value than the default 1e-6 of matlab ode solvers.
        %
        % See the AbsTol description in the odeset documentation for the effects of this property.
        %
        % @propclass{critical} Absolute error tolerances require adoption to the current situation
        % and scale of computations.
        %
        % @default 1e-9;
        %
        % See also: RelTol ode15i odeset
        AbsTol = 1e-9;
    end
        
    methods
        function this = MLode15i
            this = this@solvers.ode.BaseImplSolver;
            this.Name = 'MatLab ode15i implicit solver wrapper';
            
            this.registerProps('RelTol','AbsTol');
        end
        
        function [t,x] = implicit_solve(this, implfun, t, x0, x0p, opts)
            if nargin < 6
                opts = odeset;
            end
            opts = odeset(opts, 'RelTol', this.RelTol, 'AbsTol', this.AbsTol);
            [t,x] = ode15i(implfun, t, x0, x0p, opts);
            x = x';
            t = t';
        end
    end
    
end