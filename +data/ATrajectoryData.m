classdef ATrajectoryData < handle
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
% @new{0,5,dw,2011-10-20} Added the ATrajectoryData.getBoundingBox method and
% implementations in MemoryTrajectoryData and FileTrajectoryData.
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
 
    methods
        function transferFrom(this, source)
            % @todo write code that copies the trajectories from one
            % ATrajectoryData to another
            nt = source.getNumTrajectories;
            this.clearTrajectories;
            for nr=1:nt
                [x, mu, inputidx, ctime] = source.getTrajectoryNr(nr);
                this.addTrajectory(x, mu, inputidx, ctime);
            end
        end
        
        function [t, dia, m, M] = computeManifoldDiameters(this)
            % @todo implement!
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
        
        % Returns the degrees of freedom for the trajectories and parameter size
        [d, mud] = getTrajectoryDoFs(this);
    end
end

