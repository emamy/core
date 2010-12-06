classdef RandomSampler < sampling.BaseSampler
    %RANDOMSAMPLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The number of samples to take.
        Samples = 30;
    end
    
    methods
        function samples = performSampling(this, model)
            % Randomly generates input samples by choosing params and
            % time parameter by chance.
            % @ingroup s_rand
            %
            % Note: Read test_model.m notices on sampling settings.
            
            sys = model.System;
            
            % Compute rescaling factor for each Desired entry to end up with
            % prod(<ranges>) = model.sampling.samples
            rescaler = (this.Samples/prod([sys.Params(:).Desired]))^(1/sys.ParamCount);
            
            ranges = cell(sys.ParamCount,1);
            
            %% Compute param ranges
            for pidx=1:sys.ParamCount
                p = sys.Params(pidx);
                num = round(p.Desired * rescaler);
                ranges{pidx} = sort(rand(1,num) * (p.MaxVal-p.MinVal) + p.MinVal);
            end
            
            % No combinations necessary if only one parameter
            if sys.ParamCount == 1
                samples = ranges{pidx};
            else
                samples = general.Utils.createCombinations(ranges);
            end
        end
    end
    
end

