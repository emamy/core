classdef ModelData < data.FileData
% Data class that contains a model's large data, including subspace matrices, trajectories and
% approximation training data.
%
%
% @new{0,6,dw,2012-07-10} Created this general class with all detailed data and a root folder
% where the model's (possible) data resides. The former ATrajectoryData has been split up into this
% and ATrajectoryData.
%
% @author Daniel Wirtz @date 2012-07-10
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
        
        % The trajectory training data for the full model
        %
        % @type data.ATrajectoryData @default data.MemoryTrajectoryData
        TrajectoryData = [];
        
        % Evaluations of the system's nonlinearity on the trajectoriy data.
        %
        % @todo write setter
        %
        % @type data.FileDataCollection @default []
        TrajectoryFxiData = [];
        
        % Training data for the core function approximation.
        %
        % @type data.ApproxTrainData @default []
        ApproxTrainData = [];
        
        % Training data for the jacobian approximation.
        %
        % So far only used by the DEIM error estimator but placed here as
        % ATrajectoryData is the container for all large model offline data.
        %
        % @type data.ApproxTrainData @default []
        JacobianTrainData = [];
                
        % The projection matrix `V` for the reduced subspace.
        %
        % @type matrix<double> @default []
        V = [];
        
        % The V-biorthogonal matrix `W` for the reduced subspace (`W^tV=I_d`)
        %
        % @type matrix<double> @default []
        W = [];
        
        % A struct containing precomputed values regarding the similarity
        % transform of the jacobians in the context of the
        % error.DEIMEstimator.
        %
        % @type struct @default []
        JacSimTransData = struct;
    end
    
    properties(SetAccess=private)
        SimCache;
    end
    
    properties(Dependent)
        % The number of samples contained in the model data
        SampleCount;
    end
    
    methods
        function this = ModelData(varargin)
            % Creates a new container for large full model data.
            %
            % Parameters:
            % varargin: Either a models.BaseFullModel instance to infer the data directory
            % from, or a string containing a valid folder. @default A temporary folder within
            % the KerMor.TempDirectory
            if isempty(varargin)
                data_dir = fullfile(KerMor.App.TempDirectory,sprintf('temp_md_%s',...
                    general.IDGenerator.generateID));
            elseif isa(varargin{1},'models.BaseFullModel')
                data_dir = fullfile(KerMor.App.DataStoreDirectory,sprintf('model_%s',varargin{1}.ID));
            elseif ischar(varargin{1})
                data_dir = varargin{1};
            else
                error('Invalid argument: %s',class(varargin{1}));
            end
            this = this@data.FileData(data_dir);
            this.TrajectoryData = data.MemoryTrajectoryData;
            this.SimCache = data.FileTrajectoryData(fullfile(data_dir,'simcache'));
        end
        
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
        
        function delete(this)
            this.SimCache = [];
            this.TrajectoryData = [];
            this.ApproxTrainData = [];
            this.JacobianTrainData = [];
            delete@data.FileData(this);
        end
    end
    
    %% Getter & Setter
    methods
        function value = get.SampleCount(this)
            value = size(this.ParamSamples,2);
        end
             
        function set.ApproxTrainData(this, value)
            if ~isempty(value) && ~isa(value, 'data.ApproxTrainData')
                error('The ApproxTrainData must be a data.ApproxTrainData subclass.');
            end
            if ~isequal(this.ApproxTrainData, value)
                this.ApproxTrainData = [];
                this.ApproxTrainData = value;
            end
        end
        
        function set.TrajectoryData(this, value)
            if ~isempty(value) && ~isa(value, 'data.ATrajectoryData')
                error('The TrajectoryData must be a data.ATrajectoryData subclass.');
            end
            % Precaution: Unset the old traj data as it may trigger deletion of the same folder
            % that is used by the new one.
            if ~isequal(this.TrajectoryData,value)
                this.TrajectoryData = [];
                this.TrajectoryData = value;
            end
        end
        
        function set.JacobianTrainData(this, value)
            if ~isempty(value) && ~isa(value, 'data.ApproxTrainData')
                error('The JacobianTrainData must be a data.ApproxTrainData subclass.');
            end
            if ~isequal(this.JacobianTrainData, value)
                this.JacobianTrainData = [];
                this.JacobianTrainData = value;
            end
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

