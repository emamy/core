classdef IAutoConfig < handle
    % Interface for approximation algorithms that can "guess" their
    % configuration.
    %
    % Introduced to enable highmost automatization of the offline
    % phase. So far (at Version 0.3) the component-wise kernel
    % approximations make use of this in order to guess Gaussian kernel
    % (gamma) parameters.
    %
    % @author Daniel Wirtz @date 2011-03-21
    %
    % See also: BaseCompWiseKernelApprox
    %
    % @new{0,3,dw,2011-03-21} Added this class to KerMor. Now all
    % implementing subclasses are autoconfigured (if enabled) during the
    % offline generation process (Used in
    % models.BaseFullModel.off4_genApproximationTrainData)
    
    properties
        % Flag that indicates if automatic configuration is enabled.
        %
        % Set to false if you want to manually configure the approximation.
        % @default = true;
        EnableAutoconf = true;
    end
    
    methods
        function autoConfig(this, model)
            % Automatically detects suitable settings using the given
            % ModelData.
            %
            % Performs validity checks and calls an internal template
            % method that has to be implemented in any subclass.
            %
            % Parameters:
            % model: The Model instance associated with the current
            % model. Its Data property (Instance of models.ModelData) must
            % already have set the fields ApproxTrainData and
            % ApproxfValues.
            %
            % See also: ModelData
            % BaseFullModel.off4_genApproximationTrainData
            if ~isa(model, 'models.BaseFullModel')
                error('The model parameter has to be a models.BaseFullModel subclass.');
            end
            if isempty(model.Data.ApproxTrainData) || isempty(model.Data.ApproxfValues)
               error('The fields ApproxTrainData and ApproxfValues of the current model data have to be set.'); 
            end
            
            if this.EnableAutoconf
                this.autoconfigure(model);
            end
        end
    end
    
    methods(Access=protected, Abstract)
        % Template method for internal configuration call
        %
        % Parameters:
        % model: The Model instance. The fields ApproxTrainData and
        % ApproxfValues of the Data field are guaranteed to be set.
        autoconfigure(this, model);
    end
    
end

