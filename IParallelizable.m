classdef IParallelizable < handle
    % IPARALLELIZABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Flag whether the code should be run in parallel or not.
        ComputeParallel = false; 
    end    
    
    methods
        function IParallelizable = IParallelizable % class constructor
           t = matlabpool('size');
            if (t > 0)
                IParallelizable.ComputeParallel = true;
            else 
                disp('Proceeding with Serial Processing...');
            end
        end
    end
     
end