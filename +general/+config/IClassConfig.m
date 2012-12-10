classdef IClassConfig < handle
% IClassConfig: 
%
% @docupdate
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
        n = getNumConfigurations(this);
        
        applyConfiguration(this, nr, object);
        
        str = getConfigurationString(this, nr);
    end
    
end