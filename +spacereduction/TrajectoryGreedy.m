classdef TrajectoryGreedy < spacereduction.BaseSpaceReducer
% TrajectoryGreedy: Greedy subspace computation over a fixed set of trajectories.
%
% This subspace computation method uses a Greedy algorithm to determine a subspace which has a
% projection error less than the specified tolerance Eps.
% Therein, subsequently the largest trajectory projection error norm sums are taken for all
% trajectories, and the one with the largest is chosen. Of that trajectory, a POD of the projection
% error (orthogonal complement) vectors is performed and the first mode is added to the base, which
% is orthonormalized afterwards.
%
% @author Daniel Wirtz @date 2011-08-04
%
% @change{0,5,dw,2011-09-29} Fixed initial subspace setup when all vectors
% are identical or even zero. Also added a new alglimit-property
% TrajectoryGreedy.MaxSubspaceSize to have an upper bound to the subspace
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
        
        % The maximum subspace size
        %
        % @propclass{alglimit}
        %
        % @default 100 @type integer
        MaxSubspaceSize = 100;
    end
    
    methods
        function this = TrajectoryGreedy
            this.registerProps('Eps');
        end
        
        function [V,W] = generateReducedSpace(this, model)
            pod = general.POD;
            pod.Value = 1; % MUST stay at 1 or getInitialSpace will fail.
            pod.Mode = 'abs';
            
            o = general.Orthonormalizer;
            o.Algorithm = 'gs';
            
            md = model.Data;
            V = this.getInitialSpace(md, pod);
            [err, idx] = this.getMaxErr(V, md);
            cnt = 1; 
            % Maximum possible subspace size
            ss = size(V,1);
            while (err(end) > this.Eps) && cnt < this.MaxSubspaceSize && size(V,2) < ss
                x = md.getTrajectoryNr(idx);
                e = x - V*(V'*x);
                Vn = pod.computePOD(e);
                V = o.orthonormalize([V Vn]);
                [err(end+1), idx] = this.getMaxErr(V,md);%#ok
                cnt = cnt+1;
            end
            if cnt == this.MaxSubspaceSize
                fprintf('POD-Greedy: Maximum predefined subspace size %d reached. Aborting.\n',this.MaxSubspaceSize);
            end
            if size(V,2) == ss
                fprintf('POD-Greedy: Maximum possible subspace size %d reached. Aborting.\n',ss);
            end
            if KerMor.App.Verbose > 1
                fprintf('Trajectory greedy finished with subspace size %d and error %e.\n',size(V,2),err(end));
                if KerMor.App.Verbose > 2
                    figure;
                    semilogy(err);
                    xlabel('POD-Greedy Iterations'); ylabel('Error norms');
                    title('L_2 norm sum over whole projection error trajectory vectors');
                end
            end
            W = V;
        end
    end
    
    methods(Access=private)
        
        function V = getInitialSpace(this, md, pod)
            % Simplest: POD the first trajectory
            %x = md.getTrajectoryNr(1);
            %V = pod.computePOD(x);
            
            % More advanced: Compute first POD mode of the initial values!
            n = md.getNumTrajectories();
            x = md.getTrajectoryNr(1);
            %x0 = zeros(size(x,1),n);
            x0 = x(:,1);
            for idx=2:n
                x = md.getTrajectoryNr(idx);
                x = x(:,1);
                if isempty(general.Utils.findVecInMatrix(x0,x))
                    x0 = [x0 x];%#ok
                end
            end
            if size(x0,2) > 1
                V = pod.computePOD(x0);
            elseif all(x0 == 0)
                if KerMor.App.Verbose > 1
                    fprintf('Initial values are all zero vectors. Using main POD mode of first trajectory as initial space.\n');
                end
                V = pod.computePOD(md.getTrajectoryNr(1));
            else
                V = x0;
            end
            V = V / norm(V);
        end
        
        function [maxerr, midx] = getMaxErr(this, V, md)%#ok
            midx = -1;
            maxerr = 0;
            for k=1:md.getNumTrajectories
                x = md.getTrajectoryNr(k);
                e = x - V*(V'*x);
                e = sum(sqrt(sum(e.^2,1)));
                if maxerr < e
                    maxerr = e;
                    midx = k;
                end
            end
             if KerMor.App.Verbose > 0
                fprintf('POD-Greedy max error for subspace size %d: %e\n',size(V,2),e);
            end
        end
    end
    
    methods(Static)
        function test_TrajectoryGreedy
            
            sv = 15;
            dim = 100;
            dt = .05;
            T = 1;
            
            m = models.BaseFullModel;
            m.dt = dt;
            m.T = T;
            
            m.Sampler = sampling.GridSampler;
            
            m.Approx = [];
            
            tg = spacereduction.TrajectoryGreedy;
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