classdef ModelData < handle
    %MODELDATA Class that contains a model's large data
    %   @todo setter for matrices V,W with checks for norm one etc
    
    properties
        % A Model's parameter samples
        ParamSamples = [];
        
        % A Model's Snapshots
        Snapshots = [];
        
        % A Model's f (nonlinearity) - values at the snapshot points
        fValues = [];
        
        % The projection matrix for the reduced subspace.
        % All columns must have norm one.
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
        PlainSnapshotArray;
    end
    
    methods
        function idx = getSampleIndex(this, mu)
            % Finds the column index of the given parameter vector `\mu`
            % within the Data's ParamSamples matrix. Returns [] if `\mu` is
            % not found.
            % 
            % See also: ModelData/getTrajectory
            idx = [];
            for n=1:this.SampleCount
                if isequal(this.ParamSamples(:,n),mu)
                    idx = n;
                    return;
                end
            end
%             if ~isempty(ps)
%                 idx = find(sum(repmat(mu,1,this.SampleCount) == ps,1) == size(ps,1));
%             end
        end
        
        function x = getTrajectory(this, mu, inputidx)
            % Gets a system's trajectory for the given `\mu` and
            % inputindex.
            % Returns [] if no trajectory is found in the Data's Snapshots.
            % 
            % See also: ModelData/getSampleIndex
            x = [];
            if nargin == 2 || isempty(inputidx)
                inputidx = 1;
            end
            % Ensure that data for the given inputidx is available
            if size(this.Snapshots,4) >= inputidx
                pidx = this.getSampleIndex(mu);
                x = this.Snapshots(:,:,pidx,inputidx);
            end
        end
    end
    
    %% Getter & Setter
    methods
        function set.ParamSamples(this, value)
            % @todo Getter & setter / validation
            this.ParamSamples = value;
        end
        
        function set.Snapshots(this, value)
            this.Snapshots = value;
        end
        
        function set.fValues(this, value)
            this.fValues = value;
        end
        
%         function set.V(this, value)
%             if ~all(round(value' * value) == diag(ones(size(value,2),1)))
%                 warning('KerMor:ModelData:Vortho','V matrix may not be orthogonal!');
%             end
%             this.V = value;
%         end
        
        function value = get.SampleCount(this)
            value = size(this.ParamSamples,2);
        end
        
        function value = get.PlainSnapshotArray(this)
            value = this.Snapshots(:,:);
        end
        
    end
    
%     methods(Static)
%         function res = test_getTrajectory
%             
%         end
%         
%         function res = test_getSampleIndex
%             
%         end
%     end
    
end

