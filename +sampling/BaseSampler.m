classdef BaseSampler < KerMorObject
    %BaseSampler Basis class for parameter sampling classes
    
    properties
        % A domain from which the produced/generated samples may come from.
        % Use in subclasses at the sampling.BaseSampler.performSampling
        % method.
        %
        % @type sampling.Domain @default []
        %
        % See also: sampling.Domain
        Domain = [];
    end
    
    methods
        function samples = generateSamples(this, model)
            sys = model.System;
            if sys.ParamCount == 0
                samples = [];
            else
                params = 1:sys.ParamCount;
                if ~isempty(model.TrainingParams)
                    params = params(model.TrainingParams);
                end
                samples = this.performSampling(sys.Params(params));
                if ~isempty(model.TrainingParams)
                    full = repmat(model.DefaultMu,1,size(samples,2));
                    full(model.TrainingParams,:) = samples;
                    samples = full;
                end
            end
        end
        
        function set.Domain(this, value)
            if ~isempty(value) && ~isa(value,'sampling.Domain')
                error('The Domain property must be a sampling.Domain instance.');
            end
            this.Domain = value;
        end
    end
    
    methods(Abstract)
        % Template method for actual sampling.
        samples = performSampling(params)
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if ~isa(obj,'sampling.BaseSampler')
                if nargin < 2
                    error('Invalid call to loadobject method. Must pass the struct to instantiate the class from.');
                end
                if isfield(from,'Domain')
                    obj.Domain = from.Domain;
                else
                    obj.Domain = [];
                end
                obj = loadobj@KerMorObject(obj, from);
            else
                obj = loadobj@KerMorObject(obj);
            end
        end
    end
    
    methods(Static)
        function res = test_SubsetSampling
            % setup parameter domain etc
            % domain are all points in unit square with norm > 0.7
            m = models.BaseFullModel;
            s = models.BaseFirstOrderSystem(m);
            m.System=s;
            s.addParam('param_a',.5);
            s.addParam('param_b',1);
            s.addParam('param_c',0,'Range',[0,1]);
            s.addParam('param_d',2);
            
            m.TrainingParams = [1 3];
            m.DefaultMu = rand(4,1);
            
            m.Sampler = sampling.RandomSampler;
            m.Sampler.Samples = 500;
            samples = m.Sampler.generateSamples(m);
            
            res = size(samples,1) == 4;
            res = res & all(samples(2,:) == m.DefaultMu(2));
            res = res & all(samples(4,:) == m.DefaultMu(4));
            
            m.Sampler = sampling.GridSampler;
            samples = m.Sampler.generateSamples(m);
            res = res & size(samples,1) == 4;
            res = res & all(samples(2,:) == m.DefaultMu(2));
            res = res & all(samples(4,:) == m.DefaultMu(4));
        end
    end
end

