classdef InterpolConfig < IClassConfig
% InterpolConfig: 
%
% @author Daniel Wirtz @date 2013-01-23
%
% @new{0,7,dw,2013-01-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        
        function this = InterpolConfig
            this.RequiredPrototypeClass = 'general.interpolation.KernelInterpol';
            this.Prototype = general.interpolation.KernelInterpol;
        end
        
        function n = getNumConfigurations(~)
            n = 1;
        end
        
        
        function object = configureInstance(this, ~)
            % Returns the number of configurations that can be applied
            %
            % Return values:
            % object: The class object for which to apply the configuration @type handle
            object = this.getProtoClass;
        end
        
        function str = getConfigurationString(~, ~, ~)
            % Returns the number of configurations that can be applied
            %
            % Return values:
            % str: The configuration string @type char
            str = '';
        end
        
        function str = getConfiguredPropertiesString(~)
            str = 'none';
        end
        
        function conf = getSubPart(this, ~, ~)
            conf = this;
        end
        
    end
    
    methods(Access=protected)    
        
        function collectRanges(~,~,~)
            % do nothing
        end
    end
    
end