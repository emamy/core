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
    
    properties
       FileMatrixCacheBlockSizeGB = 1; 
    end
    
    properties(Access=private)
        model;
    end
    
    methods
        function this = VectorialSemiImplicitEuler(model)
            this.model = model;
            this.Name = 'Vectorial semi-implicit Euler''s method';
        end
    end
    
    methods(Access=protected)
        function td = customSolve(this, ~, t, x0, outputtimes)
            % Implements the actual semi-implicit solving of the given ODE system.
            %
            %
            % @change{0,7,dw,2013-01-11} Using the outputtimes parameter in order to provide a
            % more memory-efficient implementation.
            s = this.model.System;
            if isempty(s.A)
                error('This solver requires an (affine) linear system component A.');
            end
            
            % Initialize result
            steps = length(t);
            n = size(x0,1);
            starttime = tic;
            
            dt = t(2)-t(1);
            if ~isempty(find((abs(t(2:end)-t(1:end-1) - dt)) / dt > 1e-6,1)) %any(t(2:end)-t(1:end-1) - dt > 100*eps)
                error('non-equidistant dt timesteps.');
            end
            
            effsteps = length(find(outputtimes));
            numSolves = size(s.mu,2);  % number of columns of mu
            numTimestepsPerBlock = min(effsteps,floor(this.FileMatrixCacheBlockSizeGB*1024^3/(8*n*numSolves)));
            blocksize = 8*n*numSolves*numTimestepsPerBlock/(1024^2);
            x = data.FileMatrix(n,effsteps*numSolves,'BlockSize',blocksize);

            % Assign initial values
            x(:,1:numSolves) = x0;

            % Initialize output index counter
            outidx = 2;
            
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
            M = speye(n-edim);
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
                l = zeros(n,n*numSolves);
                u = zeros(n,n*numSolves);
                for index = 1:numSolves
                    [l(:,(index-1)*n+1:index*n),u(:,(index-1)*n+1:index*n)] = lu(M - dt * A(:,(index-1)*n+1:index*n));
                end
            end
            
            % Solve for each time step
            oldx = x0(1:end-edim,:);
            newx = zeros(size(oldx));
            if ~isempty(s.u)
                [urow,ucol] = size(s.u);
                if (ucol ~= numSolves) && (ucol ~= 1)
                    error('u and mu must have same number of columns or u must be column vector')
                end
            end
            for idx = 2:steps;
                RHS = M*oldx;
                if ~isempty(s.f)
                    RHS = RHS + dt*s.f.evaluate(oldx, t(idx-1), s.mu);
                end
                if ~isempty(s.u)
                    if numSolves == 1
                        RHS = RHS + dt*s.B.evaluate(t(idx-1), s.mu)*s.u(t(idx-1));
                    else
                        B = s.B.evaluate(t(idx-1), s.mu);
                        input = zeros(n,numSolves);
                        if ucol == 1   % u is column-vector
                            for index = 1:numSolves
                                input(:,index) = B(:,(index-1)*urow+1:index*urow)*s.u(t(idx-1));
                            end
                        else  % u is matrix with mucol columns
                            for index = 1:numSolves
                                input(:,index) = B(:,(index-1)*urow+1:index*urow)*s.u(t(idx-1),:,index);   % todo: indexing of u ???
                            end
                        end
                        RHS = RHS + dt*input;
                    end
                    
                end
                
                % Time-independent case: constant
                if ~fdep && ~mdep
                    for index = 1:numSolves
                        newx(:,index) = u(:,(index-1)*n+1:index*n)\(l(:,(index-1)*n+1:index*n)\RHS(:,index));
                    end
                    % Time-dependent case: compute time-dependent components with updated time
                else
                    if fdep
                        A = s.A.evaluate(1, t(idx), s.mu);
                    end
                    if mdep
                        M = s.M.evaluate(t(idx));
                    end
                    for index = 1:numSolves
                        [l,u] = lu(M - dt * A(:,(index-1)*n+1:index*n));
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
                    x(:,(outidx-1)*numSolves+1:outidx*numSolves) = [newx; newex];

                    % x(:,outidx,:) = [newx; newex];%#ok
                    % squeeze(x(:,t,:)) == x at timeindex t
                    outidx = outidx+1;
                end
                oldx = newx;
                oldex = newex;
            end
            
            % "Guess" computation time by dividing by the number of simultaneous solves
            ctime = toc(starttime)/numSolves;
            td = this.transformFileMatrixCache(x, numTimestepsPerBlock, ctime);
        end
    end
    
    methods(Access=private)
        function td = transformFileMatrixCache(this, x, numTimestepsPerBlock, ctime)
            s = this.model.System;
            numSolves = size(s.mu,2);
            
            % Extend input indices if necessary
            inidx = s.inputidx;
            if isscalar(inidx)
                inidx = repmat(inidx,1,numSolves);
            end
            muhash = general.Utils.getHash([s.mu; inidx]);
            savedir = fullfile(this.model.Data.DataDirectory,sprintf('simres_%s',muhash));
            td = data.FileTrajectoryData(savedir);
            td.UniformTrajectories = false;
            td.ReplaceExisting = true;
            
            % This "equation" holds as the FileMatrix used for caching the results was designed
            % that way.
            nCols = numSolves*numTimestepsPerBlock;
            
            for k = 1:x.getNumBlocks % Loop through block number for correct positions
                for idx = 1:numSolves % Extract timesteps for each solve
                    mu = s.mu(:,idx);
                    % These indices effectively only access one block
                    xpos = (k-1)*nCols + (idx:numSolves:nCols);
                    if k == 1 % At first block, just add as there is not existing data
                        td.addTrajectory(x(:,xpos),mu,inidx(idx),ctime);
                    else % Else: Load existing data, extend, save
                        % Caution! The last block of the FileMatrix cache may not have full
                        % size, thus cut all indices that would exceed it.
                        if k == x.getNumBlocks
                            xpos(xpos > x.m) = [];
                        end
                        % Get existing data
                        xold = td.getTrajectory(mu,inidx(idx));
                        % Extend
                        xold = [xold x(:,xpos)]; %#ok
                        % Store
                        td.addTrajectory(xold,mu,inidx(idx),ctime);
                    end
                end
            end
            
            mu_values = s.mu; input_indices = s.inputidx; themodel = this.model;%#ok
            save(savedir,'mu_values','input_indices','themodel');
            
%             Forgot that timesteps are not even stored with trajectory data!
%             t = t(outputtimes);
%             tpos = (k-1)*numTimestepsPerBlock + (1:numTimestepsPerBlock);
        end
    end
end