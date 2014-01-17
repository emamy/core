classdef IParallelizable < handle
    % IParallelizable Interface to indicate parallel computation capability
    % of the implementing classes
    %
    % @author Syed Ammar @date 2011-04-14
    %
    % @change{0,7,dw,2014-01-17} Removed the
    % KerMor.UseMatlabParallelComputing flag and hence this class is only
    % of minimal functionality now. However, having an interface to
    % indicate parallel computation capability is not making anything else
    % more complicated, which is why this remains in KerMor.
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
            this.ComputeParallel = value;
        end
    end
     
end