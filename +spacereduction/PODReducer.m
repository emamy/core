classdef PODReducer < spacereduction.BaseSpaceReducer & general.POD & general.IReductionSummaryPlotProvider
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
    
    properties
        IncludeTrajectoryFxiData = false;
    end
    
    properties(SetAccess=private)
        % The singular values of the SVD of the used trajectories.
        %
        % @type rowvec<double> @default []
        SingularValues = [];
    end
        
    methods
        function [V,W] = generateReducedSpace(this, model)
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
            
            [V, this.SingularValues] = this.computePOD(td);
            
            % Here W=V!
            W = V;
        end
        
        function plotSummary(this, pm, context)
            if ~isempty(this.SingularValues)
                str = sprintf('%s: POD singular value decay',context);
                h = pm.nextPlot('podreducer_singvals',str,...
                    'subspace size','singular values');
                semilogy(h,this.SingularValues,'LineWidth',2);
            else
                warning('spacereduction:PODReducer',...
                    'Singular value data empty. Not providing summary.');
            end
        end
    end
    
end

