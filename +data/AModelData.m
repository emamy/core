classdef AModelData < handle
% Data class that contains a model's large data, including subspace matrices, trajectories and
% approximation training data.
%
% @author Daniel Wirtz @date 2011-08-03
%
% @new{0,5,dw,2011-08-03} 
% - Added this class.
% - Removed the field ApproxfValues and put it into a struct with ApproxTrainData.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
%
% @todo Find different representation of model data / trajectory saving! (allows sequential POD of
% trajectories for subspace comp)
    
    properties
        % A Model's parameter samples
        ParamSamples = [];
        
        % Training data for the core function approximation.
        %
        % This struct has the following fields:
        % xi: The system's state variable `x_i=x(t_i)`.
        % ti: The times `t_i`
        % mui: The parameter values `\mu_i`
        % fxi: The f evaluations `f(x_i)`. If subspace projection is used, this will equal
        % `f(Vz_i)`.
        ApproxTrainData = struct('xi',[],'ti',[],'mui',[],'fxi',[]);
                
        % The projection matrix `V` for the reduced subspace.
        V;
        
        % The V-biorthogonal matrix `W` for the reduced subspace (`W^tV=I_d`)
        W;
    end
    
    properties(Dependent)
        % The number of samples contained in the model data
        SampleCount;
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
        x = getTrajectory(this, mu, inputidx);
        
        % Gets the total number of trajectories
        n = getNumTrajectories(this);
        
        % Gets the trajectory with the number nr.
        [x, mu, inputidx] = getTrajectoryNr(this, nr);
        
        % Adds a trajectory to the ModelData instance.
        addTrajectory(this, x, mu, inputidx);
        
        % Clears all stored trajectory data.
        clearTrajectories(this);
    end
    
    %% Getter & Setter
    methods
        function value = get.SampleCount(this)
            value = size(this.ParamSamples,2);
        end
             
        function set.ApproxTrainData(this, value)
            if ~isa(value, 'struct') || ~isfield(value,'xi') || ~isfield(value,'ti')...
                    || ~isfield(value,'mui') || ~isfield(value,'fxi')
                error('The property must be struct with the fields xi,ti,mui and fxi.');
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

