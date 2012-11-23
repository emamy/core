classdef EpsSVRConfig < general.config.IClassConfig
% EpsSVRConfig: 
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
        Epsilons = [];
        Lambdas = [];
    end
    
    methods
        function this = EpsSVRConfig(values)
            if nargin == 1 && size(values,1) == 2
                this.Epsilons = values(1,:);
                this.Lambdas = values(2,:);
            end
        end
        
        function n = getNumConfigurations(this)
            n = length(this.Epsilons);
        end
        
        function applyConfiguration(this, nr, svr)
            svr.Eps = this.Epsilons(nr);
            svr.Lambda = this.Lambdas(nr);
        end
        
        function str = getConfigurationString(this, nr)
            str = [];
            if ~isempty(this.Epsilons)
                str = sprintf('Eps: %g, Lambda: %g',this.Epsilons(nr),this.Lambdas(nr));
            end
        end
    end
    
end