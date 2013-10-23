classdef BaseSampler < KerMorObject
    %BaseSampler Basis class for parameter sampling classes
    
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
        % Template method for actual sampling.
        samples = performSampling(model)
    end
end

