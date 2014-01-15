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
            % Parameters:
            % model: the full or reduced model @type models.BaseModel
            %
            % Return values:
            % samples: the randomly chosen parameters, number of rows equal
            % to number of model's parameters, number of columns as specified in property Samples
            % @type matrix<double>
            %
            % @ingroup s_rand
            sys = model.System;
            if isempty(this.Seed)
                seed = round(cputime*100);
            else
                seed = this.Seed;
            end
            r = RandStream('mt19937ar','Seed',seed);
            if isempty(this.Domain)
                factor = r.rand(sys.ParamCount,this.Samples);
                miv = repmat([sys.Params(:).MinVal]',1,this.Samples);
                mav = repmat([sys.Params(:).MaxVal]',1,this.Samples);
                islog = strcmp({sys.Params(:).Spacing},'log')';
                if any(islog)
                    miv(islog,:) = log10(miv(islog,:));
                    mav(islog,:) = log10(mav(islog,:));
                end
                samples = miv + factor.*(mav-miv);
                if any(islog)
                    samples(islog,:) = 10.^samples(islog,:);
                end
                % Create samples according to the given Domain
            else
                if size(this.Domain.Points,1) ~= sys.ParamCount
                    erros('Dimension of Domain does not match the parameter count.');
                end
                samples = zeros(sys.ParamCount,this.Samples);
                pos = 1;
                while pos < this.Samples
                    % create new samples
                    factor = r.rand(sys.ParamCount,this.Samples);
                    miv = repmat([sys.Params(:).MinVal]',1,this.Samples);
                    mav = repmat([sys.Params(:).MaxVal]',1,this.Samples);
                    islog = strcmp({sys.Params(:).Spacing},'log')';
                    if any(islog)
                        miv(islog,:) = log10(miv(islog,:));
                        mav(islog,:) = log10(mav(islog,:));
                    end
                    newsamples = miv + factor.*(mav-miv);
                    if any(islog)
                        newsamples(islog,:) = 10.^newsamples(islog,:);
                    end
                    
                    % Filter the new samples by the domain and augment
                    newsamples = this.Domain.filter(newsamples);
                    numnew = size(newsamples,2);
                    spos = pos:min(pos+numnew-1, this.Samples);
                    npos = 1:min(numnew,this.Samples-pos+1);
                    samples(:,spos) = newsamples(:,npos);
                    pos = spos(end) + 1;
                end
            end
            % Sort samples if one-dimensional
            if size(samples,1) == 1
                samples = sort(samples);
            end
        end
    end
    
    methods(Static)
        
        function res = test_DomainSampling
            % setup parameter domain etc
            % domain are all points in unit square with norm > 0.7
            m = models.BaseFullModel;
            m.Sampler.domain_sampling
            s=models.BaseDynSystem(m);
            m.System=s;
            s.addParam('param_a',[0,1],10);
            s.addParam('param_b',[0,1],10);
            x = repmat((0:0.1:1),1,11);
            y = ones(11,1)*(0:0.1:1); y=y(:)';
            sampling = [x;y];
            valid = false(1,length(x));
            for i=1:length(valid)
                if norm(sampling(:,i))> 0.7
                    valid(i)=true;
                end
            end
            m.Sampler.domain_sampling = sampling;
            m.Sampler.is_valid = valid;
            m.Sampler.Samples = 50;
            samples = m.Sampler.performSampling(m);
            res = true;
            for i=1:size(samples)
                % small domain violations are possible because point with
                % norm <0.7 may be closest to a point in domain_sampling
                % with norm > 0.7
                res = res && (norm(samples(:,i)) >= 0.65);
            end
            plot(samples(1,:),samples(2,:),'o',m.Sampler.domain_sampling(1,:),m.Sampler.domain_sampling(2,:),'*',0.7*sin(0:0.05:pi/2),0.7*cos(0:0.05:pi/2))
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj)
            if ~isa(obj,'sampling.RandomSampler')
                from = obj;
                obj = sampling.RandomSampler;
                obj.Samples = from.Samples;
                obj.Seed = from.Seed;
                obj = loadobj@sampling.BaseSampler(obj,from);
            end
        end
    end
    
end

