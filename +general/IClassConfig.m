classdef IClassConfig < handle
% IClassConfig: Abstract interface for a set of configurations that can be applied to a given
% algorithm
%
% See also: kernels.RBFConfig general.regression.EpsSVRConfig
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2012-11-22} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties
        % The index of the configuration with the best results.
        %
        % Set inside the algorithm this configuration set is used for.
        %
        % @propclass{verbose} Helps identifying the effectively used configuration.
        %
        % @type integer @default []
        vBestConfigIndex = [];
    end
    
    methods(Abstract)
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % n: The number of configurations @type integer
        n = getNumConfigurations(this);
        
        % Returns the number of configurations that can be applied
        %
        % Parameters:
        % nr: The configuration number @type integer
        % object: The class object for which to apply the configuration @type handle
        applyConfiguration(this, nr, object);
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % str:  @type integer
        str = getConfigurationString(this, nr);
    end
    
end