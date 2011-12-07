classdef LinearImplEuler < solvers.ode.BaseCustomSolver
% LinearImplEuler: 
%
%
%
% @author Daniel Wirtz @date 2011-12-06
%
% @new{0,6,dw,2011-12-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        model;
    end
    
    methods
        
        function this = LinearImplEuler(model)
            this.model = model;
        end
        
    end
    
    methods(Access=protected)
        
        function x = customSolve(this, ~, t, x0)
            s = this.model.System;
%             if ~isa(s,'dscomponents.LinearCoreFun')
%                 error('boo.');
%             end
            % Initialize result
            steps = length(t);
            dt = t(2)-t(1);
            if any(t(2:end)-t(1:end-1) - dt > 100*eps)
                error('non-equidistant dt timesteps.');
            end
            
            rtm = this.RealTimeMode;
            if rtm
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                x = [x0 zeros(size(x0,1),steps-1)];
            end
            
            null = zeros(size(x0));
            A = s.f.getStateJacobian(null,0,s.mu);
            b = s.f.evaluate(null,0,s.mu);
            I = eye(size(A));
            
            % Check if a mass matrix is present
            if ~isempty(s.M)
                M = s.M.evaluate(0); 
            else
                M = I;
            end
            %[l,u] = lu(M + dt * A);
            Ai = inv(M + dt * A);
            
            % Solve for each time step
            oldx = x0;
            for idx = 2:steps;
                RHS = M*oldx + dt*b;
                if ~isempty(s.u)
                    RHS = RHS + dt*s.B.evaluate(t, s.mu)*s.u(t);
                end
                %newx = u\(l\RHS);
                newx = Ai * RHS;
                %newx = (M + dt * A)\RHS;
                
                if rtm
                    ed.Times = t(idx);
                    ed.States = newx;
                    this.notify('StepPerformed',ed);
                else
                    x(:,idx) = newx;%#ok
                end
                oldx = newx;
            end
        end
    end
    
end