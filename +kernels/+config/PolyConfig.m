classdef PolyConfig
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
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The different polynomial kernel degrees to use
        %
        % @type rowvec<double> @default []
        Degrees = [];
    end
    
    methods
        function n = getNumConfigurations(this)
            n = length(this.Degrees);
        end
        
        function applyConfiguration(this, nr, kernel)
            kernel.Degree = this.Degress(nr);
        end
        
        function str = getConfigurationString(this, nr)
            str = [];
            if ~isempty(this.Degrees)
                str = sprintf('Degree: %g',this.Degrees(nr));
            end
        end
    end
    
end