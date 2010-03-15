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
    
    methods(Access=protected)
        function comb = createCombinations(this, ranges)
            % protected function that takes an array of
            % value ranges and computes their combinations
            
            n = length(ranges);
            % Create nd-grids
            [matrices{1:n}] = ndgrid(ranges{:});
            % Convert to np x params matrix
            comb = zeros(n,numel(matrices{1}));
            for idx=1:n
                comb(idx,:) = matrices{idx}(:);
            end
        end
    end
    
end

