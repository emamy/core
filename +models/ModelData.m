classdef ModelData < handle
    % Data class that contains a model's large data
    %
    % @todo Implmement all setters with appropriate checks
    %
    % @change{0,1,dw} More common projection via matrices `V,W` instead of
    % `V,V^t`.
    % @change{0,3,dw,2011-04-01} 
    % - Changed the old 'ProjTrainData' to 'TrainingData', as this property
    % name describes the usage more precisely.
    % - Changes
    
    properties
        % A Model's parameter samples
        ParamSamples = [];
        
        % Training data for subspace & approximation computations.
        %
        % This data is the base for subspace computation algorithms located
        % in package spacereduction.
        % Each column represents a snapshot vector composed of
        % - Dim1: The index of the parameter in ParamSamples used to
        % compute this snapshot. Zero if no parameters are used.
        % - Dim2: The index of the input function used to compute this
        % snapshot. Zero if no inputs are used.
        % - Dim3: The time `t` that corresponds to the snapshot's time.
        % - Dim4-end: The system's state variable for the snapshot.
        %
        % @todo overhaul data storage format! (may be getting too big, and
        % currently replicates lots of values in the first two dimensions,
        % i.e. param and inputindices) Same for ApproxTrainData
        TrainingData = [];
        
        % Training data for the core function approximation.
        % Each vector is composed of
        % - Dim1: The index of the parameter in ParamSamples used to
        % compute this vector. Zero if no parameters are used.
        % - Dim2: The index of the input function used to compute this
        % vector. Zero if no inputs are used.
        % - Dim3: The time `t` that corresponds to `x(t)`.
        % - Dim4-end: The system's state variable for the vector, namely
        % `x(t)`.
        ApproxTrainData = [];
        
        % A Model's f (nonlinearity) - values at the snapshot points
        ApproxfValues = [];
        
        % The projection matrix for the reduced subspace.
        V;
        
        % The V-biorthogonal matrix for the reduced subspace (`W^tV=I_d`)
        W;
    end
    
    properties(Dependent)
        % The number of samples contained in the model data
        SampleCount;
        
        % Large Snapshot array.
        %
        % Returns all snapshots in an d x n matrix where d is the dimension
        % of each snapshot and n the total number of snapshots. This
        % abstracts from possibly associated parameters or inputs used for
        % generation.
        %PlainSnapshotArray;
    end
    
    methods
        
        function mu = getParams(this, idx)
            % Returns the parameter `mu` for the given indices idx. Returns
            % [] if any index is invalid or idx==[].
            mu = [];
            if ~isempty(idx) && all(idx > 0) && all(idx < size(this.ParamSamples,2)+1)
                mu = this.ParamSamples(:,idx);
            else
                %warning('models:ModelData','No parameter found for idx %d',idx);
            end
        end
        
%         function idx = getSampleIndex(this, mu)
%             % Finds the column index of the given parameter vector `\mu`
%             % within the Data's ParamSamples matrix. Returns [] if `\mu` is
%             % not found.
%             % 
%             % See also: ModelData/getTrajectory
%             idx = [];
%             for n=1:this.SampleCount
%                 if all(abs(this.ParamSamples(:,n) - mu) < sqrt(eps))
%                     idx = n;
%                     return;
%                 end
%             end
%         end
        
        % Not longer required as trajectories aren't stored completely.
%         function x = getTrajectory(this, mu, inputidx)
%             % Gets a system's trajectory for the given `\mu` and
%             % inputindex.
%             % Returns [] if no trajectory is found in the Data's Snapshots.
%             % 
%             % See also: ModelData/getSampleIndex
%             x = [];
%             if nargin == 2 || isempty(inputidx)
%                 inputidx = 1;
%             end
%             % Ensure that data for the given inputidx is available
%             if size(this.Snapshots,4) >= inputidx
%                 pidx = this.getSampleIndex(mu);
%                 x = this.Snapshots(:,:,pidx,inputidx);
%             end
%         end
    end
    
    %% Getter & Setter
    methods
        function value = get.SampleCount(this)
            value = size(this.ParamSamples,2);
        end     
    end
end

