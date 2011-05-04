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
        
    methods
        function [V,W] = generateReducedSpace(this, model)
            % Implements the abstract method from BaseSpaceReducer
            
            data = model.Data.TrainingData;   
            
            % Perform POD on state variable part of the snapshots!
            V = this.computePOD(data(4:end,:));
            
            % Here W=V!
            W = V;
        end
    end
    
end

