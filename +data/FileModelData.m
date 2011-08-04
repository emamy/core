classdef FileModelData < data.AModelData
% FileModelData: 
%
%
%
% @author Daniel Wirtz @date 2011-08-04
%
% @new{0,5,dw,2011-08-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
    end
    
    methods
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
        
        function n = getNumTrajectories(this)
           % Gets the total number of trajectories
        end
        
        function [x, mu, inputidx] = getTrajectoryNr(this, nr)
            % Gets the trajectory with the number nr.
        end
        
        function addTrajectory(this, x, mu, inputidx)
            % Adds a trajectory to the ModelData instance.
        end
        
        function clear(this)
        end
    end
    
end