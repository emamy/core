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
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(Constant, Access=private)
        FileTrajectoryDataFolder = 'trajectories';
        FileTrajectoryIncompleteDataFolder = 'trajectories_incomplete';
        FileTrajectoryFxiDataFolder = 'trajectories_fx';
        SimCacheDataFolder = 'simcache';
    end
    
    properties
        % A Model's parameter samples as column vector collection
        %
        % @type matrix @default []
        ParamSamples = [];
        
        % The trajectory training data for the full model (only complete
        % trajectories!)
        %
        % @type data.ATrajectoryData @default data.MemoryTrajectoryData
        TrajectoryData = [];
        
        % The incomplete trajectory training data for the full model is
        % saved separately
        %
        % @type data.ATrajectoryData @default data.MemoryTrajectoryData
        TrajectoryIncompleteData = [];
        
        % Evaluations of the system's nonlinearity on the trajectoriy data.
        %
        % @todo write setter
        %
        % @type data.FileDataCollection @default []
        TrajectoryFxiData = [];
        
        % The orthonormalized vectors spanning the space of the input components.
        %
        % Set the models.BaseFullModel.ComputeBSpan to true to compute during
        % models.BaseFullModel.off2_genTrainingData.
        %
        % @type matrix<double> @default []
        InputSpaceSpan = [];
        
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
        
        % Array of ProjectionSpace instances that contains all the separate
        % linear projection spaces computed during the offline phase.
        %
        % For "normal" settings, this is just one global space. For
        % multifield simulations (e.g. Biomechanics / FEM) separate
        % projection spaces for displacement and velocity could be needed.
        %
        % @type colvec<data.ProjectionSpace> @default []
        ProjectionSpaces;
        
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
            % from, or a string containing a valid folder. If a model instance is passed, and
            % additional varargin argument may contain a target directory. In this case, the
            % model's models.BaseFullModel.SaveTag along with the ID will be used to determine
            % the model data storage directory.@default A temporary folder within
            % the KerMor.TempDirectory
            %
            % @change{0,7,dw,2013-04-15} Using the SaveTag + ID to determine the data storage
            % location if a BaseFullModel is passed. Then, a second parameter may be the base
            % directory for the model's data directory.
            if isempty(varargin)
                data_dir = fullfile(KerMor.App.TempDirectory,sprintf('temp_md_%s',...
                    IDGenerator.generateID));
            elseif isa(varargin{1},'models.BaseFullModel')
                m = varargin{1};
                basedir = KerMor.App.DataDirectory;
                if length(varargin) > 1 && ischar(varargin{2})
                    basedir = varargin{2};
                end
                data_dir = fullfile(basedir,sprintf('%s_id-%s',m.SaveTag,m.ID));
            elseif ischar(varargin{1})
                data_dir = varargin{1};
            else
                error('Invalid argument: %s',class(varargin{1}));
            end
            this = this@data.FileData(data_dir);
            % Init with default memory data storages
            this.TrajectoryData = data.MemoryTrajectoryData;
            this.TrajectoryIncompleteData = data.MemoryTrajectoryData;
            this.TrajectoryFxiData = data.MemoryTrajectoryData;
            
            %% Sim cache is FileTrajectoryData by default
            this.SimCache = data.FileTrajectoryData(fullfile(data_dir,this.SimCacheDataFolder));
            % Dont force uniform trajectory lengths for simulation cache
            this.SimCache.UniformTrajectories = false;
        end
        
        function addProjectionSpace(this, V, W, dims)
            % Adds a new projection space to the model data instance.
            if ~isa(V,'data.FileMatrix')
                V = data.FileMatrix(V,'Dir',this.DataDirectory);
            end
            if ~isempty(W) && ~isa(W,'data.FileMatrix')
                W = data.FileMatrix(W,'Dir',this.DataDirectory);
            end
            s = data.ProjectionSpace(V,W,dims);
            if isempty(this.ProjectionSpaces)
                this.ProjectionSpaces = s;
            else
                this.ProjectionSpaces(end+1) = s;
            end
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
            this.ProjectionSpaces = [];
            this.SimCache = [];
            this.TrajectoryData = [];
            this.TrajectoryIncompleteData = [];
            this.TrajectoryFxiData = [];
            this.ApproxTrainData = [];
            this.JacobianTrainData = [];
            delete@data.FileData(this);
        end
        
        function relocate(this, new_dir)
            % Relocates this ModelData instances filesystem-bound quantities to a new location.
            %
            % Parameters:
            % new_dir: The new directory @type char
            %
            % @new{0,7,dw,2013-05-28} Added this method.
            if nargin < 2
                new_dir = Utils.getDir('Please specify the directory',KerMor.App.DataDirectory);
                if new_dir == 0
                    fprintf(2,'No directory selected. Aborting.\n');
                    return;
                end
            end
            
            if ~isempty(this.TrajectoryData) && isa(this.TrajectoryData,'data.FileData')
                this.TrajectoryData.relocate(...
                    fullfile(new_dir,this.FileTrajectoryDataFolder));
            end
            if ~isempty(this.TrajectoryIncompleteData) && isa(this.TrajectoryIncompleteData,'data.FileData')
                this.TrajectoryIncompleteData.relocate(...
                    fullfile(new_dir,this.FileTrajectoryIncompleteDataFolder));
            end
            if ~isempty(this.TrajectoryFxiData) && isa(this.TrajectoryFxiData,'data.FileData')
                this.TrajectoryFxiData.relocate(...
                    fullfile(new_dir,this.FileTrajectoryFxiDataFolder));
            end
            if ~isempty(this.SimCache) && isa(this.SimCache,'data.FileData')
                this.SimCache.relocate(...
                    fullfile(new_dir,this.SimCacheDataFolder));
            end
            if ~isempty(this.ApproxTrainData)
                this.ApproxTrainData.relocate(new_dir);
            end
            if ~isempty(this.JacobianTrainData)
                this.JacobianTrainData.relocate(new_dir);
            end
            if ~isempty(this.ProjectionSpaces)
                for k=1:length(this.ProjectionSpaces)
                    this.ProjectionSpaces(k).relocate(new_dir);
                end
            end
            relocate@data.FileData(this, new_dir);
        end
        
        function useFileTrajectoryData(this, overwrite)
            % Sets the TrajectoryData and TrajectoryFxiData classes to filesystem based
            % versions.
            %
            % Throws an error if any affected instance is already a data.FileTrajectoryData.
            %
            % Parameters:
            % overwrite: Set to true if new instances should be created, overwriting previous
            % FileTrajectoryData instances. @type logical @default false
            if nargin < 2
                overwrite = false;
            end
            
            % complete trajectories
            if ~overwrite && isa(this.TrajectoryData,'data.FileTrajectoryData')
                error('TrajectoryData already is a data.FileTrajectoryData. Use "true" parameter to overwrite.');
            end
            old = this.TrajectoryData;
            this.TrajectoryData = [];
            this.TrajectoryData = data.FileTrajectoryData(...
                fullfile(this.DataDirectory,this.FileTrajectoryDataFolder));
            % Copy existing data
            this.TrajectoryData.transferFrom(old);
            
            % incomplete trajectories
            if ~overwrite && isa(this.TrajectoryIncompleteData,'data.FileTrajectoryData')
                error('TrajectoryIncompleteData already is a data.FileTrajectoryData. Use "true" parameter to overwrite.');
            end
            old = this.TrajectoryIncompleteData;
            this.TrajectoryIncompleteData = [];
            this.TrajectoryIncompleteData = data.FileTrajectoryData(...
                fullfile(this.DataDirectory,this.FileTrajectoryIncompleteDataFolder));
            % Copy existing data
            this.TrajectoryIncompleteData.transferFrom(old);
            
            % trajectory fxi data
            if ~overwrite && isa(this.TrajectoryFxiData,'data.FileTrajectoryData')
                error('TrajectoryFxiData already is a data.FileTrajectoryData. Use "true" parameter to overwrite.');
            end
            old = this.TrajectoryFxiData;
            this.TrajectoryFxiData = [];
            this.TrajectoryFxiData = data.FileTrajectoryData(...
                fullfile(this.DataDirectory,this.FileTrajectoryFxiDataFolder));
            % Copy existing data
            this.TrajectoryFxiData.transferFrom(old);            
        end
        
        function [t, y, mu, inputidx, ct] = getCachedTrajectory(this, nr)
            % Shorthand for retrieval of cached trajectories that have been
            % created during simulations of a models.BaseFullModel.
            %
            % See also: models.BaseFullModel.computeTrajectory
            [y, mu, inputidx, ct] = this.SimCache.getTrajectoryNr(nr);
            t = 0:mu(2):mu(1);
            mu = mu(3:end);
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            % Loads the ModelData with automatic relocation for different hosts.
            %
            % If the directory does not exist and we are on a different host, open a directory
            % selection dialog to select the new data directory.
            if ~isa(this,'data.ModelData')
                initfrom = this;
                this = data.ModelData(initfrom.fDataDir);
                this = loadobj@data.FileData(this,initfrom);
                if isfield(initfrom,'V') && ~isempty(initfrom.V)
                    this.addProjectionSpace(initfrom.V,initfrom.W,size(initfrom.V,1));
                end
            else
                this = loadobj@data.FileData(this);
            end
            if exist(this.DataDirectory,'file') == 0
                fprintf(2,'If you have moved the folder or attempt to load data on a different machine, use the ModelData.relocate method.\n');
            end
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
        
        function set.ProjectionSpaces(this, value)
            if ~isempty(value) && ~isa(value,'data.ProjectionSpace')
                error('Value has to be a data.ProjectionSpace instance/array');
            end
            this.ProjectionSpaces = value;
        end
        
    end
end

