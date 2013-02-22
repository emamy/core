classdef RandomSampler < sampling.BaseSampler
    %RandomSampler Selects Samples many random parameters.
    %
    % @author Daniel Wirtz @date 2011-10-11
    %
    % @new{0,7,dw,2013-02-22} Added a new property Spacing to also allow for logarithmical
    % sampling over the parameter ranges.
    %
    % @new{0,5,dw,2011-10-11} Added this class (former version renamed to
    % WeightedRandomSampler)
    
    properties(SetObservable)
        % The number of samples to take.
        %
        % @propclass{critical} Determines how many parameter samples are taken and thus directly
        % the offline computation time and model approximation quality.
        %
        % @default 30
        Samples = 30;
        
        % The seed for the random number generator.
        %
        % If empty, a seed computed from the cputime is used.
        % 
        % @propclass{optional} Set to empty for new choices on every run.
        %
        % @type integer @default []
        Seed = [];
        
        % The spacing method.
        %
        % Available options are 'log' and 'lin' for logarithmic/linear grid
        % sampling.
        %
        % @propclass{important} Different model/parameters require
        % differently scaled sampling
        %
        % @type char @default lin
        Spacing = 'lin';
    end    
    
    methods
        
        function this = RandomSampler
            this.registerProps('Samples');
        end
        
        function set.Samples(this, value)
            if ~isposintscalar(value)
                error('Value must be a positive integer scalar');
            end
            this.Samples = value;
        end
        
        function samples = performSampling(this, model)
            % Randomly generates input samples by choosing params and
            % time parameter by chance.
            %
            % @ingroup s_rand
            sys = model.System;
            if isempty(this.Seed)
                seed = round(cputime*100);
            else
                seed = this.Seed;
            end
            r = RandStream('mt19937ar','Seed',seed);
            factor = r.rand(sys.ParamCount,this.Samples);
            miv = repmat([sys.Params(:).MinVal]',1,this.Samples);
            mav = repmat([sys.Params(:).MaxVal]',1,this.Samples);
            if strcmp(this.Spacing,'log')
                miv = log10(miv);
                mav = log10(mav);
            end
            samples = miv + factor.*(mav-miv);
            if strcmp(this.Spacing,'log')
                samples = 10.^samples;
            end
            % Sort samples if one-dimensional
            if size(samples,1) == 1
                samples = sort(samples);
            end
        end
    end
    
end

