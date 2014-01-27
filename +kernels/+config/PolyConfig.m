classdef PolyConfig < IClassConfig
% PolyConfig: Configuration settings for polynomial kernels
%
% @docupdate
%
% @author Daniel Wirtz @date 2012-11-26
%
% @new{0,7,dw,2012-11-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % The different polynomial kernel degrees to use
        %
        % @type rowvec<double> @default []
        Degrees = [];
    end
    
    methods
        
        function this = PolyConfig(values)
            if nargin == 1
                this.Degrees = values;
            end
            this.RequiredPrototypeClass = 'kernels.PolyKernel';
            this.Prototype = kernels.PolyKernel;
        end
        
        function n = getNumConfigurations(this)
            n = length(this.Degrees);
        end
        
        function k = configureInstance(this, nr)
            k = this.getProtoClass;
            k.Degree = this.Degrees(nr);
        end
        
        function str = getConfigurationString(this, nr, ~)
            str = [];
            if ~isempty(this.Degrees)
                str = sprintf('Degree: %g',this.Degrees(nr));
            end
        end
        
        function str = getConfiguredPropertiesString(~)
            str = 'Degree';
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            idx = this.getPartIndices(partNr, totalParts);
            conf = kernels.config.PolyConfig(this.Degrees(idx));
        end
        
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            this.addRange(ptable, [proppath {'Degrees'}],min(this.Degrees),max(this.Degrees));
        end
    end
    
end