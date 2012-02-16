classdef AModelData < handle
% Data class that contains a model's large data, including subspace matrices, trajectories and
% approximation training data.
%
% @author Daniel Wirtz @date 2011-08-03
%
% @change{0,6,dw,2011-12-14} Now also storing the computation time \c ctime for a trajectory.
% This comes from the fact that for error estimator comparisons we need the \c ctime for each
% trajectory for every different error estimator. See the changes in models.BaseModel for more
% information.
%
% @new{0,5,dw,2011-10-20} Added the AModelData.getBoundingBox method and
% implementations in MemoryModelData and FileModelData.
%
% @new{0,5,dw,2011-08-03} 
% - Added this class.
% - Removed the field ApproxfValues and put it into a struct with
% ApproxTrainData.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties
        % A Model's parameter samples as column vector collection
        %
        % @type matrix @default []
        ParamSamples = [];
        
        % Training data for the core function approximation.
        %
        % @type data.ApproxTrainData @default []
        ApproxTrainData = [];
                
        % The projection matrix `V` for the reduced subspace.
        %
        % @type matrix @default []
        V = [];
        
        % The V-biorthogonal matrix `W` for the reduced subspace (`W^tV=I_d`)
        %
        % @type matrix @default []
        W = [];
    end
    
    properties(Dependent)
        % The number of samples contained in the model data
        SampleCount;
    end
    
    methods
        function mu = getParams(this, idx)
            % Returns parameters for given indices.
            %
            % Parameters:
            % idx: The parameter indices @type rowvec<int32>
            %
            % Return values: 
            % mu: The parameters `mu` for the given indices idx. Returns
            % [] if any index is invalid or idx==[]. @type matrix<double>
            mu = [];
            if ~isempty(idx) && all(idx > 0) && all(idx < size(this.ParamSamples,2)+1)
                mu = this.ParamSamples(:,idx);
            else
                %warning('models:ModelData','No parameter found for idx %d',idx);
            end
        end
        
        function idx = getSampleIndex(this, mu)
            % Finds the column index of the given parameter vector `\mu`
            % within the Data's ParamSamples matrix. Returns [] if `\mu` is
            % not found.
            % 
            % See also: ModelData/getTrajectory
            idx = [];
            for n=1:this.SampleCount
                if all(abs(this.ParamSamples(:,n) - mu) < sqrt(eps))
                    idx = n;
                    return;
                end
            end
        end
    end
    
    methods(Abstract)
        % Gets the traejctory for the given parameter `\mu` and input index.
        [x, ctime] = getTrajectory(this, mu, inputidx);
        
        % Gets the total number of trajectories
        n = getNumTrajectories(this);
        
        % Gets the trajectory with the number nr.
        [x, mu, inputidx, ctime] = getTrajectoryNr(this, nr);
        
        % Adds a trajectory to the ModelData instance.
        addTrajectory(this, x, mu, inputidx, ctime);
        
        % Clears all stored trajectory data.
        clearTrajectories(this);
        
        % Gets the bounding box of the state space of all trajectories.
        [x,X] = getBoundingBox(this);
    end
    
    %% Getter & Setter
    methods
        function value = get.SampleCount(this)
            value = size(this.ParamSamples,2);
        end
             
        function set.ApproxTrainData(this, value)
            if ~isempty(value) && ~isa(value, 'data.ApproxTrainData')
                error('The property must be a data.ApproxTrainData subclass.');
            end
            this.ApproxTrainData = value;
        end
        
        function set.V(this, value)
            if ~isa(value, 'double')
                error('value must be a valid matrix of type double');
            end
            this.V = value;
        end
        
        function set.W(this, value)
            if ~isa(value, 'double')
                error('value must be a valid matrix of type double');
            end
            this.W = value;
        end
        
    end
    
end

