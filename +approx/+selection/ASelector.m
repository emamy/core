classdef ASelector < KerMorObject & ICloneable
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
    % @new{0,4,dw,2011-05-04} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @new{0,3,dw,2011-04-12} Added this class to implement strategy
    % pattern for training data selection.
    
    properties(SetAccess=protected)
        % The indices of the selected approximation training data of the last call to
        % selectTrainingData.
        %
        % @propclass{data}
        %
        % @default []
        LastUsed = [];
    end
    
    methods
        
        function this = ASelector
            this = this@KerMorObject;
            this.registerProps('LastUsed');
        end
        
        function copy = clone(this, copy)
            copy.LastUsed = this.LastUsed;
        end
        
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
            
            if isempty(this.LastUsed)
                warning('KerMor:TrainDataSelection','TrainDataSelector did not set the LastUsed property.\nComputing automatically, which ist very slow for large data sets!');
                this.LastUsed = general.Utils.findVecInMatrix(model.Data.TrainingData,atd);
                %error('The subclass does not set the LastUsed property after selection. Please implement in order for other KerMor classes to work properly.');
            end
        end
        
        function set.LastUsed(this, value)
            if ~isposintmat(value)
                error('The LastUsed property must contain a matrix with positive integer indices');
            end
            this.LastUsed = value;
        end
    end   
    
    methods(Abstract, Access=protected)
        % Template method for subclasses to specify selection behaviour.
        atd = select(this, model)
    end
end

