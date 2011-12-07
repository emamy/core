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
        
        function x = customSolve(this, ~, t, x0)
            s = this.model.System;
%             if ~isa(s,'dscomponents.LinearCoreFun')
%                 error('boo.');
%             end
            % Initialize result
            steps = length(t);
            dt = t(2:end)-t(1:end-1);
            
            rtm = this.RealTimeMode;
            if rtm
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                x = [x0 zeros(size(x0,1),steps-1)];
            end
            
            A = s.f.evaluate(0,x0);
            I = eye(size(A));
            % Solve for each time step
            for idx = 2:steps;
                
                LHS = M + dt(idx-1) * A;
                % Check if a mass matrix is present
                if ~isempty(s.M)
                    M = s.M.evaluate(t); 
                else
                    M = I;
                end
                RHS = M*x(:,idx-1);
                if ~isempty(s.u)
                    RHS = RHS + s.B.evaluate(t, s.mu)*s.u(t);
                end
                newx = LHS\RHS;
                
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