classdef ScalarNuSVRConfig < general.IClassConfig
    % ScalarNuSVRConfig:
    %
    % @docupdate
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
    
    properties(SetAccess=private)
        nuvals = [];
    end
    
    methods
        function this = ScalarNuSVRConfig(nuvals)
            if nargin == 1
                this.nuvals = nuvals;
            end
        end
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % n: The number of configurations @type integer
        function n = getNumConfigurations(this)
            n = length(this.nuvals);
        end
        
        % Returns the number of configurations that can be applied
        %
        % Parameters:
        % nr: The configuration number @type integer
        % object: The class object for which to apply the configuration @type handle
        function applyConfiguration(this, nr, object)
            object.nu = this.nuvals(nr);
        end
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % str:  @type integer
        function str = getConfigurationString(this, nr)
            str = sprintf('nu=%g',this.nuvals(nr));
        end
        
    end
    
    methods(Access=protected)        
        function collectRanges(this, proppath)
            this.addRange([proppath {'Nu'}],min(this.nuvals),max(this.nuvals));
        end        
    end
    
end