classdef WeightedRandomSampler < sampling.RandomSampler
    %WeightedRandomSampler: Computes random samples using the Desired
    %fields of the parameters.
    %
    % @docupdate
    %
    % @author Daniel Wirtz @date 2010-10-11
    %
    % @new{0,5,dw,2011-10-11} Renamed the former RandomSampler to this
    % class as it is not truly random.
    %
    % @change{0,3,sa,2011-05-10} Implemented property setter
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    
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
                samples = Utils.createCombinations(ranges);
            end
        end
    end
    
end

