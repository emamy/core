classdef PODReducer < spacereduction.BaseSpaceReducer & general.POD
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
    % @change{0,5,dw,2011-08-04} Adopted to the new ModelData structure. Now the global POD array is
    % assembled from all the trajectories available.
        
    methods
        function [V,W] = generateReducedSpace(this, model)
            % Implements the abstract method from BaseSpaceReducer
            
            % Compile all trajectories into a large vector
            nt = model.Data.getNumTrajectories;
            if nt == 0
                error('No training data found in ModelData.');
            end
            
            x = model.Data.getTrajectoryNr(1);
            tlen = size(x,2);
            data = zeros(size(x,1),tlen*nt);
            for k=1:nt
                data(:,(k-1)*tlen+1:k*tlen) = model.Data.getTrajectoryNr(k);
            end
            
            % Perform POD on state variable part of the snapshots!
            V = this.computePOD(data);
            
            % Here W=V!
            W = V;
        end
    end
    
end

