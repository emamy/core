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
    % @change{0,5,dw,2011-11-09} Changed the return type of selectTrainingData to an instance
    % of data.ApproxTrainData (not having the fxi property set yet)
    %
    % @new{0,5,dw,2011-08-04} Removed the LastUsed property as it is incompatible with the new
    % data.AModelData structure (the approximation training data is not a subset of TrainingData
    % anymore). Subclasses have been changed in order to adopt to the new structure, too.
    %
    % @new{0,4,dw,2011-05-04} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @new{0,3,dw,2011-04-12} Added this class to implement strategy
    % pattern for training data selection.
    %
    % @todo move to +data/+selection package
    
    methods
        
        function this = ASelector
            this = this@KerMorObject;
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
            % atd: The approximation training data @type data.ApproxTrainData
            if ~isa(model,'models.BaseFullModel')
                error('The model parameter must be a BaseFullModel subclass.');
            elseif model.Data.getNumTrajectories == 0
                error('No training data available to select approximation training data from.');
            end
            
            [xi, ti, mui] = this.select(model);
            atd = data.ApproxTrainData(xi, ti, mui);
        end
    end
    
    methods(Abstract, Access=protected)
        % Template method for subclasses to specify selection behaviour.
        %
        % Parameters:
        % model: The full model with the training data @type models.BaseFullModel
        %
        % Return values:
        % xi: The selected `x_i = x(t_i)` training data @type matrix
        % ti: The selected training times `t_i` @type rowvec
        % mui: The selected parameter samples `\mu_i` with which the states
        % `x_i` have been reached @type matrix
        [xi, ti, mui] = select(this, model)
    end
end

