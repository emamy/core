classdef PODReducer < spacereduction.BaseSpaceReducer & general.POD & IReductionSummaryPlotProvider
    %PODREDUCER Uses POD for reduced space generation.
    %
    % Internally the SVD decomposition of the snapshot array is used.
    % Several modes are supported to enable more specific reduced space
    % selection.
    %
    % See also: general.POD.Mode general.POD.Value
    %
    % @author Daniel Wirtz @date 19.03.2010
    %
    % @change{0,6,dw,2012-07-13} Now also works with data.ABlockedData arguments as trajectory
    % data.
    %
    % @change{0,5,dw,2011-08-04} Adopted to the new ModelData structure.
    % Now the global POD array is assembled from all the trajectories available.
    
    properties(SetAccess=private)
        % The singular values of the SVD of the used trajectories.
        %
        % @type rowvec<double> @default []
        SingularValues = [];
    end
    
    properties
        % Flag indicating if to include the initial space 
        % 
        % type logical @type false
        IncludeInitialSpace = false;
    end
        
    methods
        function [V, W] = generateReducedSpace(this, model)
            % Implements the abstract method from BaseSpaceReducer
            
            % Compile all trajectories into a large vector
            md = model.Data;
            td = md.TrajectoryData;
            nt = td.getNumTrajectories;
            if nt == 0
                error('No training data found in ModelData.');
            end
            
            % Augment block data with fxi values
            if this.IncludeTrajectoryFxiData
                if isempty(md.TrajectoryFxiData)
                    error('No training fxi data found in ModelData.');
                end
                td = data.JoinedBlockData(td, md.TrajectoryFxiData);
            end
            % Wrap in finite difference adder
            if this.IncludeFiniteDifferences
                td = data.FinDiffBlockData(td);
            end
            Vex = [];
            if this.IncludeBSpan
                Vex = md.InputSpaceSpan;
            end
            if this.IncludeInitialSpace
                is = this.getInitialSpace(md.TrajectoryData);
                Vex = [Vex is];
            end
            
            [V, this.SingularValues] = this.computePOD(td, Vex);
            this.ProjectionError = flipud(cumsum(flipud(this.SingularValues)));
            if ~isempty(Vex)
                o = general.Orthonormalizer;
                V = o.orthonormalize([Vex V]);
            end
            
            % Galerkin projection
            W = [];
        end
        
        function plotSummary(this, pm, context)
            plotSummary@spacereduction.BaseSpaceReducer(this, pm, context);
            if ~isempty(this.SingularValues)
                str = sprintf('%s: POD singular value decay',context);
                h = pm.nextPlot('podreducer_singvals',str,...
                    'subspace size','singular values');
                semilogy(h,this.SingularValues,'LineWidth',2);
%                 h = pm.nextPlot('podreducer_projerr',str,...
%                     'subspace size','projection error');
%                 semilogy(h,this.SingularValues,'LineWidth',2);
            else
                warning('spacereduction:PODReducer',...
                    'Singular value data empty. Not providing summary.');
            end
        end
    end
    
    methods(Access=private)
        
        function V = getInitialSpace(this, md)
            % Computes the initial space, which is the first POD mode of
            % the initial values!
            
            pod = general.POD;
            pod.Value = 1; % MUST stay at 1 or getInitialSpace will fail.
            pod.Mode = 'abs';
            n = md.getNumBlocks;
            x = md.getBlock(1);
            x0 = x(:,1);
            for idx=2:n
                    x = md.getBlock(idx);
                    x = x(:,1);
                    % Only add nonexisting vectors
                    if isempty(Utils.findVecInMatrix(x0,x))
                        x0 = [x0 x];%#ok
                    end
            end
            if all(x0(:) == 0)
                if KerMor.App.Verbose > 1
                    fprintf('Initial values are all zero vectors. Using main POD mode of first block data as initial space.\n');
                end
                V = pod.computePOD(md.getBlock(1));
            elseif size(x0,2) > 1
                V = pod.computePOD(x0);
            else
                V = x0;
            end
            V = V / norm(V);
        end
    end
    
end

