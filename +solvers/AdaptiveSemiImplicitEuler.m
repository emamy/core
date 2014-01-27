classdef AdaptiveSemiImplicitEuler < solvers.BaseCustomSolver
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
    
    properties
        epsilon;    % (exp(L(b-a))-1)/L * epsilon (Interval (a,b), Lipschiz constant L) is upper error bound
        q;
        h0;
    end
    
    methods
        function this = AdaptiveSemiImplicitEuler(model, epsilon, q, h0)
            
            if nargin < 4
                h0 = 1e-3;
                if nargin < 3
                    q = 0.3;
                    if nargin < 2
                        epsilon = 1e-2;
                    end
                end
            end
            this.h0 = h0;
            this.q = q;
            this.epsilon = epsilon;
            this.model = model;
            this.Name = 'adaptive semi-implicit euler method';
        end
    end
    
    methods(Access=protected)
        function x = customSolve(this, ~, t, x0, ~)
            % Implements the actual semi-implicit solving of the given ODE system.
            %
            % Parameters:
            % t: Either a two dimensional vector with t(1) < t(2) specifiying
            % the start and end time, or a greater or equal to three
            % dimensional, strictly monotoneously increasing vector explicitly
            % setting the desired output times. Depending on the MaxStep
            % property, the solver can work with a finer time step internally.
            % x0: The initial value
            %
            % Return values:
            % x: The computed trajectory @type matrix<double>
            
            s = this.model.System;
            if isempty(s.A)
                error('This solver requires an (affine) linear system component A.');
            end
            
            h = this.h0;
            time = t(1);
            t_end = t(end);
            
            
            rtm = this.RealTimeMode;
            
            steps = 0;
            % Create return matrix in size of effectively desired timesteps
            x = zeros(size(x0,1),length(t));
            x(:,1) = x0;
            oldx = x0;
            % Initialize output index counter
            outidx = 2;
            
            oldex = []; newex = []; edim = 0; est = [];
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
                [l1, u1] = lu(M - h * A);
            end
            
            if ~mdep && ~fdep
                [l2, u2] = lu(M - h/2 * A);
            end
            
            % Solve for each time step
            %oldx = x0(1:end-edim);
            count = 0;
            while time < t_end;
                delta = 0.5 * this.epsilon * (1-this.q);   % garantees that the inner while loop is entered
                firstloop = true;
                %time
                
                while ((delta/this.epsilon) < (1 - this.q)) || ((delta/this.epsilon) > (1 + this.q))
                    if ~firstloop
                        h = h * (this.epsilon / delta);
                        count = count +1;
                        %h
                        %delta 
                        %outidx
                        %time
                    else
                        firstloop = false;
                    end
                    x1 = oldx;
                    x2 = oldx;
                    
                    steps = steps +1;
                    
                    RHS1 = M*x1;
                    if ~isempty(s.f)
                        RHS1 = RHS1 + h*s.f.evaluate(x1, time, s.mu);
                    end
                    if ~isempty(s.u)
                        RHS1 = RHS1 + h*s.B.evaluate(time, s.mu)*s.u(time);
                    end
                    
                    % Time-independent case: constant
                    if ~fdep && ~mdep
                        x1 = u1\(l1\RHS1);
                        % Time-dependent case: compute time-dependent components with updated time
                    else
                        if fdep
                            A = s.A.evaluate(1, time, s.mu);
                        end
                        if mdep
                            M = s.M.evaluate(time);
                        end
                        x1 = (M - h * A)\RHS1;
                    end
                    
                    % the same for 2 steps of half size
                    for i=1:2
                        RHS2 = M*x2;
                        if ~isempty(s.f)
                            RHS2 = RHS2 + h/2*s.f.evaluate(x2, time, s.mu);
                        end
                        if ~isempty(s.u)
                            RHS2 = RHS2 + h/2*s.B.evaluate(time, s.mu)*s.u(time);
                        end
                        
                        % Time-independent case: constant
                        if ~fdep && ~mdep
                            x2 = u2\(l2\RHS2);
                            % Time-dependent case: compute time-dependent components with updated time
                        else
                            if fdep
                                A = s.A.evaluate(1, time, s.mu);
                            end
                            if mdep
                                M = s.M.evaluate(time);
                            end
                            x2 = (M - h/2 * A)\RHS2;
                        end
                    end
                    
                    delta = norm((x1-x2)./x2)/(2*h);                    
                    %delta = norm(x1-x2)/(h);
                end
                time = time + h;
                oldx = x2;
                % Implicit error estimation computation  - not implemented
                %                 if ~isempty(oldex)
                %                     ut = [];
                %                     if ~isempty(s.u)
                %                         ut = s.u(t(idx));
                %                     end
                %                     al = est.getAlpha(newx, t(idx), s.mu, ut);
                %                     bet = est.getBeta(newx, t(idx), s.mu);
                %
                %                     % Explicit
                %                     %newex = oldex + dt * est.evalODEPart([newx; oldex], t(idx-1), s.mu, ut);
                %                     newex = oldex + dt*(bet*oldex + al);
                %
                %                     % Implicit
                %                     %newex = (dt*al+oldex)/(1-dt*bet);
                %
                %                     %                     fun = @(y)y-dt*est.evalODEPart([oldx; y], t(idx), s.mu, ut)-oldex/dt;
                %                     %                     opts = optimset('Display','off');
                %                     %                     newex2 = fsolve(fun, oldex, opts);
                %                 end
                
                % Only produce output at wanted timesteps
                while time < t(end) && time > t(outidx)
                    x(:,outidx) = x2;
                    outidx = outidx + 1;
                end
            end
            l = length(t) - outidx;
            x(:,outidx:end) = repmat(x2,1,l+1);
            steps
            count
        end
    end
end