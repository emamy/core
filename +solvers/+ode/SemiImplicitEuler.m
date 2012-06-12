classdef SemiImplicitEuler < solvers.ode.BaseCustomSolver
    % SemiImplicitEuler: Solves ODEs in KerMor using implicit euler for the
    % linear part and explicit euler for the nonlinear part.
    %
    % The existence of a nonlinear part for this solver is optional. The
    % only requirement is to have the system component models.BaseDynSys.A
    %
    % If an error estimator is set for this solver, the error estimation
    % auxiliary ODE will be solved using an explicit euler scheme on the
    % way.
    %
    % This implementation replaces the former LinearImplEuler solver.
    %
    % @author Daniel Wirtz @date 2012-05-24
    %
    % @change{0,6,dw,2012-05-28} Added support for DEIMErrorEstimator use
    % with this solver.
    %
    % @new{0,6,dw,2012-05-24} Added this class.
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
        function this = SemiImplicitEuler(model)
            this.model = model;
            this.Name = 'Semi-implicit euler method';
        end
    end
    
    methods(Access=protected)
        function x = customSolve(this, ~, t, x0)
            s = this.model.System;
            if isempty(s.A)
               error('This solver requires an (affine) linear system component A.'); 
            end
%             if isa(s.Model,'models.ReducedModel') && ~isempty(s.Model.ErrorEstimator)...
%                     && s.Model.ErrorEstimator.Enabled && s.Model.ErrorEstimator.ExtraODEDims > 0
%                 error('Cannot use this solver with reduced models that have an error estimator enabled with ExtraODEDims > 0 (this solver overrides the odefun handle)');
%             end
            % Initialize result
            steps = length(t);
            dt = t(2)-t(1);
            if ~isempty(find((abs(t(2:end)-t(1:end-1) - dt)) / dt > 1e-6,1)) %any(t(2:end)-t(1:end-1) - dt > 100*eps)
                error('non-equidistant dt timesteps.');
            end
            
            rtm = this.RealTimeMode;
            if rtm
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                x = [x0 zeros(size(x0,1),steps-1)];
            end
            
            oldex = []; newex = []; edim = 0;
            if isa(this.model,'models.ReducedModel')
                est = this.model.ErrorEstimator;
                if ~isempty(est) && est.Enabled
                    edim = est.ExtraODEDims;
                    oldex = x0(end-edim+1:end);
                end
            end
            
            % Check if a mass matrix is present, otherwise assume identity matrix
            M = speye(length(x0)-edim);
            mdep = false;
            if ~isempty(s.M)
                mdep = s.M.TimeDependent;
                if ~mdep
                    M = s.M.evaluate(0);
                end
            end
            
            fdep = s.A.TimeDependent;
            if ~fdep
                % Evaluation with x=1 "extracts" the matrix A of the (affine) linear system
                A = s.A.evaluate(1, 0, s.mu);
            end
            
            % Precompute lu decomp for iteration steps if all components are not time-dependent
            if ~mdep && ~fdep
                [l, u] = lu(M - dt * A);
            end
            
            % Solve for each time step
            oldx = x0(1:end-edim);
            for idx = 2:steps;
                RHS = M*oldx;
                if ~isempty(s.f)
                    RHS = RHS + dt*s.f.evaluate(oldx, t(idx-1), s.mu);
                end
                if ~isempty(s.u)
                    RHS = RHS + dt*s.B.evaluate(t(idx-1), s.mu)*s.u(t(idx-1));
                end
                
                % Time-independent case: constant
                if ~fdep && ~mdep
                    newx = u\(l\RHS);
                    % Time-dependent case: compute time-dependent components with updated time
                else
                    if fdep
                        A = s.A.evaluate(1, t(idx-1), s.mu);
                    end
                    if mdep
                        M = s.M.evaluate(t(idx-1));
                    end
                    newx = (M - dt * A)\RHS;
                end
                
                % Explicit error estimation computation
                if ~isempty(oldex)
                    ut = [];
                    if ~isempty(s.u)
                        ut = s.u(t(idx));
                    end
                    % Explicit
                    %newex = oldex + dt * est.evalODEPart([oldx; oldex], t(idx-1), s.mu, ut);
                    % Implicit
                    al = est.getAlpha([oldx; oldex], t(idx), s.mu, ut);
                    bet = est.getBeta([oldx; oldex], t(idx), s.mu);
                    newex = (dt*al+oldex)...
                        /(1-dt*bet);
%                     fun = @(y)y-dt*est.evalODEPart([oldx; y], t(idx), s.mu, ut)-oldex/dt;
%                     opts = optimset('Display','off');
%                     newex = fsolve(fun, oldex, opts);
                end
                
                % Real time mode: Fire StepPerformed event
                if rtm
                    ed.Times = t(idx);
                    ed.States = [newx; newex];
                    this.notify('StepPerformed',ed);
                    % Normal mode: Collect solution in result vector
                else
                    x(:,idx) = [newx; newex];%#ok
                end
                oldx = newx;
                oldex = newex;
            end
        end
    end
end