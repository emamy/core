classdef MLImplSolver < solvers.BaseSolver
% MLImplSolver: 
%
%
%
% @author Daniel Wirtz @date 2011-04-14
%
% @new{0,3,dw,2011-04-14} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
        
    methods
        function this = MLImplSolver
            this.Name = 'MatLab ode15i implicit solver wrapper';
        end
        
        function [t,y] = solve(this, odefun, t, x0)
            implfun = @(t,y,yp)yp - odefun(t,y);
            yp0 = odefun(0,x0);
            
            opts = [];
            if ~isempty(this.MaxStep)
                opts = odeset(opts, 'MaxStep', this.MaxStep);
            end
            if ~isempty(this.InitialStep)
                opts = odeset(opts, 'InitialStep', this.InitialStep);
            end
            
            [t,y] = ode15i(implfun,t,x0,yp0,opts);
            
            y = y';
            t = t';
        end
    end
    
end