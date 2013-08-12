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
        
        % Samples from the full parameter domain box specified by the
        % parameter's MinVal and MaxVal . domain_sampling(:,isvalid)
        % specifies the domain in which the parameter samples must lie.
        % Dimension: `d_p \times n_s`, where `n_s` is the number of samples
        % and `d_p` is the dimension of the parameter space, i.e. the
        % system's ParamCount.
        %
        % @type matrix<double> @default []
        domain_sampling = [];
        
        % Logical vector specifying which of the points in domain_sampling
        % are valid. Must have the same number of columns as
        % domain_sampling.
        %
        % @type vector<logical> @default []
        is_valid = []
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
            if isempty(this.domain_sampling)
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
            else
                if size(this.domain_sampling,2) ~= size(this.is_valid,2)
                    error('Properties domain_sampling and is_valid must have the same number of columns.');
                elseif ~islogical(this.isvalid)
                    error('Property is_valid must be a logical vector.');
                elseif size(this.domain_sampling,1) ~= sys.ParamCount
                    erros('Number of columns of domain_sampling does not match System.ParamCount.');
                end
                %                 samples = zeros(sys.ParamCount,this.Samples);
                %                 need_resampling = 1:this.Samples;
                %                 while ~isempty(need_resampling)
                %                     % create as many new samples as needed
                %                     len = length(need_resampling);
                %                     factor = r.rand(sys.ParamCount,len);
                %                     miv = repmat([sys.Params(:).MinVal]',1,len);
                %                     mav = repmat([sys.Params(:).MaxVal]',1,len);
                %                     islog = strcmp({sys.Params(:).Spacing},'log')';
                %                     if any(islog)
                %                         miv(islog,:) = log10(miv(islog,:));
                %                         mav(islog,:) = log10(mav(islog,:));
                %                     end
                %                     newsamples = miv + factor.*(mav-miv);
                %                     if any(islog)
                %                         newsamples(islog,:) = 10.^newsamples(islog,:);
                %                     end
                %                     % check, which of the new samples are valid
                %                     idx = dsearchn(this.domain_sampling', newsamples'); % those indices of domain_sampling that are closest to new_samples
                %                     valid_idx = this.is_valid(idx);
                %                     samples(:,valid_idx)=newsamples(:,valid_idx);
                %                     need_resampling = need_resampling(~valid_idx);
                %                 end
                samples = zeros(sys.ParamCount,this.Samples);
                good_samples = 0;
                while true
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
                    % check, which of the new samples are valid
                    idx = dsearchn(this.domain_sampling', newsamples'); % those indices of domain_sampling that are closest to new_samples
                    valid_idx = this.is_valid(idx);
                    good_new_samp = sum(valid_idx);
                    if good_samples + good_new_samp < this.Samples  % need another while loop
                        samples(:,good_samples+1:good_samples+good_new_samp) = newsamples(:,valid_idx);
                        good_samples= good_samples+good_new_samp;
                    else % enough, add only as many valid samples as needed
                        newidx=1:this.Samples;
                        newidx=newidx(valid_idx);
                        newidx = newidx(1:this.Samples-good_samples);
                        samples(:,good_samples+1:end)=newsamples(:,newidx);
                        break
                    end
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
    
end

