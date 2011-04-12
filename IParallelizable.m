classdef IParallelizable < handle
    % IPARALLELIZABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Flag whether the code should be run in parallel or not.
        ComputeParallel = false; 
    end    
    
    methods
        function set.ComputeParallel(this, h)
            if h == true
                a = KerMor.App;
                if a.UseMatlabParallelComputing == true
                    this.ComputeParallel = h;
                else
                    disp('WARNING : Matlab Parallel Processing toolbox not activated');
                end
            end    
        end
    end
     
end