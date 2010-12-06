classdef GridSampler < sampling.BaseSampler
    %GRIDSAMPLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function samples = performSampling(this, model)%#ok
            % Uses the given model and generates a training set by creating
            % a regular grid in joint time/parameter space
            % @ingroup s_grid
            sys = model.System;
            ranges = cell(sys.ParamCount,1);
            
            % Create linearly spaced parameter value ranges
            for pidx=1:sys.ParamCount
                % Cater for single parameters 
                if (sys.Params(pidx).MinVal == sys.Params(pidx).MaxVal)
                    ranges{pidx} = sys.Params(pidx).MinVal;
                else
                    ranges{pidx} = linspace(sys.Params(pidx).MinVal,...
                        sys.Params(pidx).MaxVal,...
                        sys.Params(pidx).Desired);
                end
            end
            
            % Return first range if only one available.
            if sys.ParamCount == 1
                samples = ranges{1};
            else
                samples = general.Utils.createCombinations(ranges);
            end
        end
    end
    
end

