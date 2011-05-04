classdef BaseSampler < KerMorObject
    %BASESAMPLER Summary of this class goes here
    %   Detailed explanation goes here
    
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

