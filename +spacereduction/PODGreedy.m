classdef PODGreedy < spacereduction.BaseSpaceReducer & IParallelizable
% PODGreedy: Greedy subspace computation over a fixed set of trajectories.
%
% This subspace computation method uses the POD-Greedy algorithm to determine a subspace which has a
% projection error less than the specified tolerance Eps.
% Therein, subsequently the largest trajectory projection error norm sums are taken for all
% trajectories, and the one with the largest is chosen. Of that trajectory, a POD of the projection
% error (orthogonal complement) vectors is performed and the first mode is added to the base, which
% is orthonormalized afterwards.
%
% See http://www.agh.ians.uni-stuttgart.de/publications/2011/Haa11/ and
% references therein for details.
%
% @author Daniel Wirtz @date 2011-08-04
%
% @change{0,5,dw,2011-10-11} 
% - Renamed this class to PODGreedy from TrajectoryGreedy.
% - Added another stopping criterion
% PODGreedy.MinRelImprovement in order to stop subspace computations as
% soon as the subsequent improvement becomes too small.
%
% @change{0,5,dw,2011-09-29} Fixed initial subspace setup when all vectors
% are identical or even zero. Also added a new alglimit-property
% PODGreedy.MaxSubspaceSize to have an upper bound to the subspace
% size.
%
% @new{0,5,dw,2011-08-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The desired accuracy over all training trajectories
        %
        % @propclass{critical} Ultimately determines the projection error on all training
        % trajectories.
        %
        % @default 1e-6 @type double
        Eps = 1e-6;
        
        % The minimum required relative error improvement compared to the
        % previous error
        %
        % @propclass{important} Too little minimal improvement causes a too large
        % subspace
        %
        % @default 0.005 @type double
        MinRelImprovement = .005;
        
        % The maximum subspace size
        %
        % @propclass{alglimit}
        %
        % @default 100 @type integer
        MaxSubspaceSize = 100;
    end
    
    properties(SetAccess=private)
        Errors;
    end
    
    methods
        function this = PODGreedy
            this.registerProps('Eps','MinRelImprovement','MaxSubspaceSize');
        end
        
        function [V,W] = generateReducedSpace(this, model)
            md = model.Data.TrajectoryData;
            
            if KerMor.App.Verbose > 2
                fprintf('POD-Greedy: Starting subspace computation using %d trajectories...\n',md.getNumTrajectories);
            end
            
            pod = general.POD;
            pod.Value = 1; % MUST stay at 1 or getInitialSpace will fail.
            pod.Mode = 'abs';
            
            o = general.Orthonormalizer;
            o.Algorithm = 'gs';
            
            if KerMor.App.Verbose > 3
                fprintf('POD-Greedy: Computing initial space...\n');
            end
            V = this.getInitialSpace(md, pod);
            [err, idx] = this.getMaxErr(V, md);
            cnt = 1; 
            % Maximum possible subspace size
            ss = size(V,1);
            olderr = err;
            impr = 1;
            while (err(end) > this.Eps) && cnt < this.MaxSubspaceSize && size(V,2) < ss && impr(end) > this.MinRelImprovement
                x = md.getTrajectoryNr(idx);
                e = x - V*(V'*x);
                Vn = pod.computePOD(e);
                V = o.orthonormalize([V Vn]);
                [err(end+1), idx] = this.getMaxErr(V,md);%#ok
                impr(end+1) = (olderr-err(end))/olderr;%#ok
                olderr = err(end);
                cnt = cnt+1;
            end
            this.Errors = [err; impr];
            if cnt == this.MaxSubspaceSize
                fprintf('POD-Greedy: Maximum predefined subspace size %d reached. Aborting.\n',this.MaxSubspaceSize);
            end
            if size(V,2) == ss
                fprintf('POD-Greedy: Maximum possible subspace size %d reached. Aborting.\n',ss);
            end
            if impr(end) <= this.MinRelImprovement
                fprintf('POD-Greedy: Minimum relative improvement of %f not satisfied. Aborting.\n',this.MinRelImprovement);
            end
            if KerMor.App.Verbose > 1
                fprintf('POD-Greedy: Finished with subspace size %d and error %e.\n',size(V,2),err(end));
                if KerMor.App.Verbose > 2
                    figure; subplot(1,2,1);
                    semilogy(err);
                    xlabel('POD-Greedy step'); ylabel('Error norms');
                    title('L_2 norm sum over whole projection error trajectory vectors');
                    subplot(1,2,2);
                    plot(impr);
                    xlabel('POD-Greedy step'); ylabel('Relative improvement');
                    title('Error improvement relative to previous error');
                end
            end
            W = V;
        end
    end
    
    methods(Access=private)
        
        function V = getInitialSpace(this, md, pod)
            % Computes the initial space, which is the first POD mode of
            % the initial values!
            
            n = md.getNumTrajectories();
            x = md.getTrajectoryNr(1);
            x0 = x(:,1);
            if this.ComputeParallel
                parfor idx=2:n
                    x = md.getTrajectoryNr(idx);%#ok
                    x0 = [x0, x(:,1)];
                end
                x0 = unique(x0','rows')';
            else
                for idx=2:n
                    x = md.getTrajectoryNr(idx);
                    x = x(:,1);
                    % Only add nonexisting vectors
                    if isempty(general.Utils.findVecInMatrix(x0,x))
                        x0 = [x0 x];%#ok
                    end
                end
            end
            if all(x0(:) == 0)
                if KerMor.App.Verbose > 1
                    fprintf('Initial values are all zero vectors. Using main POD mode of first trajectory as initial space.\n');
                end
                V = pod.computePOD(md.getTrajectoryNr(1));
            elseif size(x0,2) > 1
                V = pod.computePOD(x0);
            else
                V = x0;
            end
            V = V / norm(V);
        end
        
        function [maxerr, midx] = getMaxErr(this, V, md)
            midx = -1;
            maxerr = 0;
            if this.ComputeParallel
                n = md.getNumTrajectories;
                if KerMor.App.Verbose > 3
                   fprintf('POD-Greedy: Computing maximum error over %d trajectories on %d workers for subspace size %d...\n',n,matlabpool('size'),size(V,2));
                end
                err = zeros(1,n);
                parfor i=1:n
                    x = md.getTrajectoryNr(i); %#ok<PFBNS>
                    hlp = sum(sqrt(sum((x - V*(V'*x)).^2,1)));
                    if KerMor.App.Verbose > 4
                       fprintf('POD-Greedy: Error for trajectory %d: %e\n',i,hlp);
                    end
                    err(i) = hlp;
                end
                [maxerr, midx] = max(err);
            else
                for k=1:md.getNumTrajectories
                    x = md.getTrajectoryNr(k);
                    e = sum(Norm.L2(x - V*(V'*x)));
                    if maxerr < e
                        maxerr = e;
                        midx = k;
                    end
                end
            end
            if KerMor.App.Verbose > 0
                fprintf('POD-Greedy: Max error with subspace size %d: %e\n',size(V,2),maxerr);
            end
        end
    end
    
    methods(Static)
        function test_PODGreedy
            
            sv = 15;
            dim = 100;
            dt = .05;
            T = 1;
            
            m = models.BaseFullModel;
            m.dt = dt;
            m.T = T;
            
            m.Sampler = sampling.GridSampler;
            
            m.Approx = [];
            
            tg = spacereduction.PODGreedy;
            tg.Eps = 1e-5;
            m.SpaceReducer = tg;
            
            m.ODESolver = solvers.ode.ExplEuler(dt);
            
            s = models.BaseDynSystem(m);
            s.B = [];
            s.addParam('mu1',[0 1],4);
            s.addParam('mu2',[0 1],5);
            s.x0 = dscomponents.ConstInitialValue(rand(dim,1)*5);
            s.MaxTimestep = dt;
            m.System = s;
            
            f = dscomponents.ParamTimeKernelCoreFun;
            f.Kernel = kernels.GaussKernel(40);
            f.TimeKernel = kernels.NoKernel;
            f.ParamKernel = kernels.GaussKernel(2);
            f.Ma = rand(dim,sv);
            f.Centers.xi = repmat(linspace(0,10,sv),dim,1);
            f.Centers.ti = 1:sv;
            f.Centers.mui = rand(2,sv);
            s.f = f;
            
            m.offlineGenerations;
        end
    end
    
end