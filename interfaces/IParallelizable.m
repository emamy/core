classdef IParallelizable < handle
    % IPARALLELIZABLE Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @author Syed Ammar @date 2011-04-14
    %
    % @change{0,3,sa,2011-04-14} Implemented UseMatlabParallelComputing functionality
    
    properties(SetObservable)
        % Flag whether the code should be run in parallel or not.
        ComputeParallel = false; 
    end    
    
    methods
        function set.ComputeParallel(this, value)
            if ~islogical(value)
                error('Value must be logical');
            end
            if value
                if KerMor.App.UseMatlabParallelComputing
                    this.ComputeParallel = value;
                else
                    warning('MatlabParallelComputing:NOT_ENABLED','Not setting flag as Matlab Parallel Processing toolbox not activated. Refer to the property KerMor.App.UseMatlabParallelComputing for more information.');
                end
            else
                this.ComputeParallel = value;
            end    
        end
    end
     
end