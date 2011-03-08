classdef ModelData < handle
    %MODELDATA Class that contains a model's large data
    %   @todo setter for matrices V,W with checks for norm one etc
    %
    % @change{0,1} More common projection via matrices `V,W` instead of
    % `V,V^t`.
    
    properties
        % A Model's parameter samples
        ParamSamples = [];
        
        % A model's snapshots for each parameter sample. This data is the
        % base for subspace computation algorithms located in package
        % spacereduction.
        % Each column represents a snapshot vector composed of
        % - Dim1: The index of the parameter in ParamSamples used to
        % compute this snapshot. Zero if no parameters are used.
        % - Dim2: The index of the input function used to compute this
        % snapshot. Zero if no inputs are used.
        % - Dim3: The time `t` that corresponds to the snapshot's time.
        % - Dim4-end: The system's state variable for the snapshot.
        ProjTrainData = [];
        
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
%         function set.ParamSamples(this, value)
%             % @todo Getter & setter / validation
%             this.ParamSamples = value;
%         end
%         
%         function set.ProjTrainData(this, value)
%             this.ProjTrainData = value;
%         end
%         
%         function set.ApproxfValues(this, value)
%             this.ApproxfValues = value;
%         end
                
        function value = get.SampleCount(this)
            value = size(this.ParamSamples,2);
        end
        
%         function value = get.PlainSnapshotArray(this)
%             value = this.Snapshots(:,:);
%         end
        
    end
    
    %         function set.V(this, value)
    %             if ~all(round(value' * value) == diag(ones(size(value,2),1)))
    %                 warning('KerMor:ModelData:Vortho','V matrix may not be orthogonal!');
    %             end
    %             this.V = value;
    %         end
    
end

