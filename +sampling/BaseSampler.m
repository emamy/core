classdef BaseSampler < handle
    %BASESAMPLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function samples = generateSamples(this, model)
            if model.System.ParamCount == 0
                samples = [];
            else
                samples = this.performSampling(model);
            end
        end
    end
    
    methods(Abstract)
        samples = performSampling(model)
    end
end

