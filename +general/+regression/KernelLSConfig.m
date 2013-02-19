classdef KernelLSConfig < general.IClassConfig
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
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetAccess=private)
        Lambdas = [];
    end
    
    methods
        function this = KernelLSConfig(lambdas)
            if nargin == 1
                this.Lambdas = lambdas;
            end
        end
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % n: The number of configurations @type integer
        function n = getNumConfigurations(this)%#ok
            n = 1;
        end
        
        % Returns the number of configurations that can be applied
        %
        % Parameters:
        % nr: The configuration number @type integer
        % object: The class object for which to apply the configuration @type handle
        function applyConfiguration(this, nr, object)%#ok
            % to nothing
        end
        
        % Returns the number of configurations that can be applied
        %
        % Return values:
        % str:  @type integer
        function str = getConfigurationString(this, nr)%#ok
            str = '';
        end
        
    end
    
    methods(Access=protected)        
        function collectRanges(this, ptable, proppath)
            this.addRange(ptable, [proppath {'Lambdas'}],min(this.Dimension),max(this.Lambdas));
        end        
    end
    
end