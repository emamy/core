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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        
        % The maximum subspace size, EXCLUDING any non-reducable dimensions
        %
        % @propclass{alglimit}
        %
        % @default 100 @type integer
        MaxSubspaceSize = 100;
        
        % An initial space to use.
        %
        % Optional, otherwise an initial space is computed using the first vectors from the
        % available trajectories.
        %
        % Must match the system's state space trajectory size.
        %
        % @propclass{optional} Can be useful if a specific choice for initial spaces is
        % available.
        %
        % @type matrix<double> @default []
        InitialSpace = [];
    end
    
    properties(SetAccess=private)
        ErrorImprovements;
    end
    
    methods
        function this = PODGreedy
            this.registerProps('Eps','MinRelImprovement','MaxSubspaceSize','InitialSpace');
        end
        
        function plotSummary(this, pm, context)
            if nargin < 2
                pm = PlotManager;
                pm.LeaveOpen = true;
                context = 'POD-Greedy';
            end
            plotSummary@spacereduction.BaseSpaceReducer(this, pm, context);
            if ~isempty(this.ErrorImprovements)
                str = sprintf('%s: Error decay over training data\nEps:%g, MaxSize:%d, FXI: %d, FD: %d, BSpan: %d',...
                    context,this.Eps,this.MaxSubspaceSize,this.IncludeTrajectoryFxiData,...
                    this.IncludeFiniteDifferences,this.IncludeBSpan);
                h = pm.nextPlot('podgreedy_gains',str,'subspace size','Gain');
                semilogy(h,this.ErrorImprovements,'LineWidth',2);
            else
                warning('spacereduction:PODGreedy',...
                    'Error data empty. Not providing summary.');
            end
        end
        
    end
    
    methods(Access=protected)
        
        function [V,W] = generateReducedSpaceImpl(this, model)
            bdata = model.Data.TrajectoryData;
            
            % Wrap in finite difference adder
            if this.IncludeFiniteDifferences
                bdata = data.JoinedBlockData(bdata, data.FinDiffBlockData(bdata));
            end
            if this.IncludeAxData
                bdata = data.JoinedBlockData(bdata, ...
                    data.AxBlockData(bdata, model.System.A, model.scaledTimes));
            end
            % Augment block data with fxi values
            if this.IncludeTrajectoryFxiData
                if isempty(model.Data.TrajectoryFxiData)
                    error('No training fxi data found in ModelData.');
                end
                bdata = data.JoinedBlockData(bdata, model.Data.TrajectoryFxiData);
            end
            
            if KerMor.App.Verbose > 2
                fprintf('POD-Greedy: Starting subspace computation using %d trajectories...\n',model.Data.TrajectoryData.getNumTrajectories);
            end
            
            reducable = this.ReducableDims;
            
            pod = general.POD;
            pod.Value = 1; % MUST stay at 1 or getInitialSpace will fail.
            pod.Mode = 'abs';
            
            o = general.Orthonormalizer;
            o.Algorithm = 'gs';
            
            % Compute initial space
            if KerMor.App.Verbose > 2
                fprintf('POD-Greedy: Computing initial space...\n');
            end
            if ~isempty(this.InitialSpace)
                V = this.InitialSpace(reducable,:);
                if this.IncludeBSpan
                    V = o.orthonormalize([V bdata.InputSpaceSpan(reducable,:)]);
                end
            elseif this.IncludeBSpan
                V = model.Data.InputSpaceSpan(reducable,:);
            else
                V = this.getInitialSpace(bdata, pod, reducable);
            end
            % Compute initial space errors
            err = [];
            for k=1:size(V,2)
               [err(end+1), idx] = this.getMaxErr(V(:,1:k), bdata, reducable);%#ok
            end
            
            cnt = 1; 
            % Maximum possible subspace size
            ss = size(V,1);
            olderr = err;
            impr = 1;
            while (err(end) > this.Eps) && cnt < this.MaxSubspaceSize && size(V,2) < ss && impr(end) > this.MinRelImprovement
                x = bdata.getBlock(idx); % get trajectory
                x = x(reducable,:);
                e = x - V*(V'*x);
                Vn = pod.computePOD(e);
                V = o.orthonormalize([V Vn]);
                [err(end+1), idx] = this.getMaxErr(V, bdata, reducable);%#ok
                impr(end+1) = (olderr-err(end))/olderr;%#ok
                olderr = err(end);
                cnt = cnt+1;
            end
            this.ErrorImprovements = impr;
            this.ProjectionError = err;
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
            % Galerkin projection!
            W = [];
        end        
        
    end
    
    methods(Access=private)
        
        function [maxerr, midx] = getMaxErr(this, V, md, reducable)
            midx = -1;
            maxerr = 0;
            if this.ComputeParallel
                n = md.getNumBlocks;
                if KerMor.App.Verbose > 3
                   fprintf('POD-Greedy: Computing maximum error over %d trajectories on %d workers for subspace size %d...\n',n,matlabpool('size'),size(V,2));
                end
                err = zeros(1,n);
                parfor i=1:n
                    x = md.getBlock(i); %#ok<PFBNS>
                    x = x(reducable,:);
                    hlp = sum(Norm.L2(x - V*(V'*x)));
                    if KerMor.App.Verbose > 4
                       fprintf('POD-Greedy: Error for block %d: %e\n',i,hlp);
                    end
                    err(i) = hlp;
                end
                [maxerr, midx] = max(err);
            else
                for k=1:md.getNumBlocks;
                    x = md.getBlock(k);
                    x = x(reducable,:);
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
            
            m.ODESolver = solvers.ExplEuler(dt);
            
            s = models.BaseDynSystem(m);
            s.B = [];
            s.addParam('mu1',.5,'Range',[0 1],'Desired',4);
            s.addParam('mu2',.6,'Range',[0 1],'Desired',5);
            s.x0 = dscomponents.ConstInitialValue(rand(dim,1)*5);
            s.MaxTimestep = dt;
            m.System = s;
            
            f = dscomponents.ParamTimeKernelCoreFun(m.System);
            kexp = kernels.ParamTimeKernelExpansion;
            kexp.Kernel = kernels.GaussKernel(40);
            kexp.TimeKernel = kernels.NoKernel;
            kexp.ParamKernel = kernels.GaussKernel(2);
            kexp.Ma = rand(dim,sv);
            kexp.Centers.xi = repmat(linspace(0,10,sv),dim,1);
            kexp.Centers.ti = 1:sv;
            kexp.Centers.mui = rand(2,sv);
            f.Expansion = kexp;
            s.f = f;
            
            m.offlineGenerations;
        end
    end
    
end