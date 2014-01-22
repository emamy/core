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
            if model.System.ParamCount == 0
                samples = [];
            else
                samples = this.performSampling(model);
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
        samples = performSampling(model)
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if nargin < 2
                error('Invalid call to loadobject method. Must pass the struct to instantiate the class from.');
            end
            if isfield(from,'Domain')
                obj.Domain = from.Domain;
            else
                obj.Domain = [];
            end
            obj = loadobj@KerMorObject(obj, from);
        end
    end
end

