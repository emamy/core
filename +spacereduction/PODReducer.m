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
        
    methods(Access=protected)        
        function [V,W] = generateReducedSpaceImpl(this, model)
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
                td = data.JoinedBlockData(td, md.TrajectoryData);
            end
            % Wrap in finite difference adder
            if this.IncludeFiniteDifferences
                td = data.FinDiffBlockData(td);
            end
            
            reducable = this.ReducableDims;
            
            Vex = [];
            if this.IncludeBSpan
                Vex = md.InputSpaceSpan(reducable,:);
            end
            if this.IncludeInitialSpace
                is = this.getInitialSpace(md.TrajectoryData, this, reducable);
                Vex = [Vex is];
            end
            
            [V, this.SingularValues] = this.computePOD(td, Vex, reducable);
            this.ProjectionError = flipud(cumsum(flipud(this.SingularValues)));
%             if ~isempty(Vex)
%                 o = general.Orthonormalizer;
%                 V = o.orthonormalize([Vex V]);
%             end
            if ~isempty(Vex)
                o = general.Orthonormalizer;
                V = [Vex V];
                id = speye(size(V,2));
                while true  % V isn't really orthogonal after orthonormalizing once
                    diff = max(max(abs(V'*V-id)));
                    if diff < 1e-9
                        break
                    else
                        V = o.orthonormalize(V);
                    end
                end
            end            
            % Galerkin projection
            W = [];
        end
    end
    
end

