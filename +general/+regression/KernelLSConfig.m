classdef KernelLSConfig < IClassConfig
% KernelLSConfig: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-01-23
%
% @new{0,7,dw,2013-01-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    properties(SetAccess=private)
        Lambdas = [];
    end
    
    methods
        function this = KernelLSConfig(lambdas)
            if nargin == 1
                this.Lambdas = lambdas;
            end
            this.RequiredPrototypeClass = 'general.regression.KernelLS';
            this.Prototype = general.regression.KernelLS;
        end
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % n: The number of configurations @type integer
        function n = getNumConfigurations(this)%#ok
            n = 1;
        end
        
        % Gets a new configured instance for a give config number
        %
        % Parameters:
        % nr: The configuration number @type integer
        %
        % Return values:
        % object: The new configured instance @type ICloneable
        function object = configureInstance(this, nr)
            object = this.getProtoClass;
            object.lambda = this.Lambdas(nr);
        end
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % str:  @type integer
        function str = getConfigurationString(this, nr)
            str = sprintf('\lambda: %g',this.Lambdas(nr));
        end
        
        function str = getConfiguredPropertiesString(~)
            str = 'lambda';
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = this.clone;
            conf.Lambdas = this.Lambdas(this.getPartIndices(partNr, totalParts));
        end
        
        function copy = clone(this)
            copy = general.regression.KernelLSConfig;
            copy = clone@IClassConfig(this, copy);
            copy.Lambdas = this.Lambdas;
        end
        
    end
    
    methods(Access=protected)        
        function collectRanges(this, ptable, proppath)
            this.addRange(ptable, [proppath {'lambda'}],min(this.Dimension),max(this.Lambdas));
        end
    end
    
end