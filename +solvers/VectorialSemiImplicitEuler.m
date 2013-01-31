classdef VectorialSemiImplicitEuler < solvers.BaseCustomSolver
    
    % not working yet!!!
    
    % vectorialSemiImplicitEuler: Solves ODEs in KerMor using implicit euler for the
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
        function this = VectorialSemiImplicitEuler(model)
            this.model = model;
            this.Name = 'vectorial Semi-implicit euler method';
        end
    end
    
    methods(Access=protected)
        function x = customSolve(this, ~, t, x0, outputtimes)
            % Implements the actual semi-implicit solving of the given ODE system.
            %
            %
            % @change{0,7,dw,2013-01-11} Using the outputtimes parameter in order to provide a
            % more memory-efficient implementation.
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
            if rtm  % not implemented
                ed = solvers.ode.SolverEventData;
                x = [];
            else
                effsteps = length(find(outputtimes));
                % Create return matrix in size of effectively desired timesteps
                % x = [x0 zeros(size(x0,1),effsteps-1)];
                % x(:,1,1) = states at t=1, mu(1),
                % x(1,:,1) = trajectory of state(1), mu(1)
                % x(1,1,:) = state(1) at t=1 for all mu
                [dim1,dim3] = size(x0);
                dim2 = effsteps;
                mucol = size(s.mu,2);  % number of columns of mu
                t_per_block = floor(0.5*1024^3 /(8*dim1/dim3));   % blocksize approx. 0.5 GB
                blocksize = 8*dim1*dim3*t_per_block/(1024^2);
                x = data.FileMatrix(dim1,dim2*dim3,'BlockSize',blocksize);
                %                 key = struct('type',{'()'},'subs',[]);
                %                 key.subs = {':', 1:dim3};
                %                 x.subsasgn(key,x0);
                x(:,1:dim3) = x0;
                % x = zeros(dim1,dim2,dim3);
                % x(:,1,:) = x0;
                % Initialize output index counter
                outidx = 2;
            end
            
            oldex = []; newex = []; edim = 0; est = [];
            if isa(this.model,'models.ReducedModel')
                est = this.model.ErrorEstimator;
                if ~isempty(est) && est.Enabled
                    edim = est.ExtraODEDims;
                    oldex = x0(end-edim+1:end,:);
                end
            end
            
            % Check if a mass matrix is present, otherwise assume identity matrix
            % M is independent from mu
            M = speye(size(x0,1)-edim);
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
                % lu decomposition in for loop in for loop over blocks of M-dt*A, significantly faster than lu decomposition of everything
                % at once
                l = zeros(dim1,dim1*mucol);
                u = zeros(dim1,dim1*mucol);
                for index = 1:mucol
                    [l(:,(index-1)*dim1+1:index*dim1),u(:,(index-1)*dim1+1:index*dim1)] = lu(M - dt * A(:,(index-1)*dim1+1:index*dim1));
                end
            end
            
            % Solve for each time step
            oldx = x0(1:end-edim,:);
            newx = zeros(size(oldx));
            if ~isempty(s.u)
                [urow,ucol] = size(s.u);
                if (ucol ~= mucol) && (ucol ~= 1)
                    error('u and mu must have same number of columns or u must be column vector')
                end
            end
            for idx = 2:steps;
                RHS = M*oldx;
                if ~isempty(s.f)
                    RHS = RHS + dt*s.f.evaluate(oldx, t(idx-1), s.mu);
                end
                if ~isempty(s.u)
                    if mucol == 1
                        RHS = RHS + dt*s.B.evaluate(t(idx-1), s.mu)*s.u(t(idx-1));
                    else
                        B = s.B.evaluate(t(idx-1), s.mu);
                        input = zeros(dim1,mucol);
                        if ucol == 1   % u is column-vector
                            for index = 1:mucol
                                input(:,index) = B(:,(index-1)*urow+1:index*urow)*s.u(t(idx-1));
                            end
                        else  % u is matrix with mucol columns
                            for index = 1:mucol
                                input(:,index) = B(:,(index-1)*urow+1:index*urow)*s.u(t(idx-1),:,index);   % todo: indexing of u ???
                            end
                        end
                        RHS = RHS + dt*input;
                    end
                    
                end
                
                % Time-independent case: constant
                if ~fdep && ~mdep
                    for index = 1:mucol
                        newx(:,index) = u(:,(index-1)*dim1+1:index*dim1)\(l(:,(index-1)*dim1+1:index*dim1)\RHS(:,index));
                    end
                    % Time-dependent case: compute time-dependent components with updated time
                else
                    if fdep
                        A = s.A.evaluate(1, t(idx), s.mu);
                    end
                    if mdep
                        M = s.M.evaluate(t(idx));
                    end
                    for index = 1:mucol
                        [l,u] = lu(M - dt * A(:,(index-1)*dim1+1:index*dim1));
                        newx(:,index) = u\(l\RHS(:,index));
                    end
                end
                
                % Implicit error estimation computation
                if ~isempty(oldex)
                    ut = [];
                    if ~isempty(s.u)
                        ut = s.u(t(idx));
                    end
                    al = est.getAlpha(newx, t(idx), s.mu, ut);
                    bet = est.getBeta(newx, t(idx), s.mu);
                    
                    % Explicit
                    %newex = oldex + dt * est.evalODEPart([newx; oldex], t(idx-1), s.mu, ut);
                    newex = oldex + dt*(bet*oldex + al);
                    
                    % Implicit
                    %newex = (dt*al+oldex)/(1-dt*bet);
                    
                    %                     fun = @(y)y-dt*est.evalODEPart([oldx; y], t(idx), s.mu, ut)-oldex/dt;
                    %                     opts = optimset('Display','off');
                    %                     newex2 = fsolve(fun, oldex, opts);
                end
                
                % Only produce output at wanted timesteps
                if outputtimes(idx)
                    if rtm  % not implemented
                        % Real time mode: Fire StepPerformed event
                        ed.Times = t(idx);
                        ed.States = [newx; newex];
                        this.notify('StepPerformed',ed);
                        % Normal mode: Collect solution in result vector
                    else
                        x(:,(outidx-1)*dim3+1:outidx*dim3) = [newx; newex];
                        
                        % x(:,outidx,:) = [newx; newex];%#ok
                        % squeeze(x(:,t,:)) == x at timeindex t
                        outidx = outidx+1;
                    end
                end
                oldx = newx;
                oldex = newex;
            end
            
            blocksize = 8*dim1*dim2/(1024^2);
            trajectory = data.FileMatrix(dim1,dim2*dim3,'BlockSize',blocksize);
            
            for idx = 1:dim3
                indices = idx:dim3:dim3*dim2;
                trajectory(:,(idx-1)*dim2+1:idx*dim2) = x(:,indices);
            end
            x = trajectory;
        end
    end
end