classdef ASelector < handle
    % Base interface for any approximation training data selection algorithm
    %
    % @attention Any subclass must set the LastUsed property to the indices corresponding to the
    % selected samples from the training data. Otherwise an error will be thrown and the generation
    % process is terminated.
    % 
    % @note Note that the selected training data is projected into the
    % precomputed subspace if spacereduction is performed.
    % 
    % @author Daniel Wirtz @date 2011-04-12
    %
    % @change{0,3,sa,2011-04-20} Implemented Setter for the property
    % LastUsed
    %
    % @new{0,3,dw,2011-04-12} Added this class to implement strategy
    % pattern for training data selection.
    
    properties(SetAccess=protected)
        % The indices of the selected approximation training data of the last call to
        % selectTrainingData.
        %
        % @default []
        LastUsed = [];
    end
    
    methods
        function atd = selectTrainingData(this, model)
            % Performs the selection procedure
            % 
            % After some validity checks the template method is called to
            % start the actual selection algorithm.
            %
            % Parameters:
            % model: The models.BaseFullModel subclass containing the
            % training data in its Data property.
            %
            % Return values:
            % atd: The approximation training data
            if ~isa(model,'models.BaseFullModel')
                error('The model parameter must be a BaseFullModel subclass.');
            elseif isempty(model.Data.TrainingData)
                error('No training data available to select approximation training data from.');
            end
            atd = this.select(model);
            
            % too slow on large sets!
            %this.LastUsed = general.Utils.findVecInMatrix(model.Data.TrainingData,atd);
            if isempty(this.LastUsed)
                error('The subclass does not set the LastUsed property after selection. Please implement in order for other KerMor classes to work properly.');
            end
        end
        
        function set.LastUsed(this, value)
            if isempy(value)
                error('The LastUsed property is empty by default. Please change in order for other KerMor classes to work properly.');
            end
            this.LastUsed = value;
        end
    end   
    
    methods(Abstract, Access=protected)
        % Template method for subclasses to specify selection behaviour.
        atd = select(this, model)
    end
end

