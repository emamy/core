classdef GaussConfig < general.config.IClassConfig
% GaussConfig: 
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
    
    properties(SetAccess=private)
        Gammas;
        Distances;
        DistEps;
    end
    
    methods
        function this = GaussConfig(varargin)
            i = inputParser;
            i.addParamValue('G',1);
            i.addParamValue('D',[]);
            i.addParamValue('Eps',eps);
            i.parse(varargin{:});
            r = i.Results;
            if ~isempty(r.D)
                ke = kernels.GaussKernel;
                g = zeros(size(r.D));
                for k = 1:length(r.D)
                    g(k) = ke.setGammaForDistance(r.D(k),r.Eps);
                end
                this.Gammas = g;
                this.Distances = r.D;
                this.DistEps = r.Eps;
            else
                this.Gammas = r.G;
                this.DistEps = [];
            end
        end
        
        function n = getNumConfigurations(this)
            if ~isempty(this.Gammas)
                n = length(this.Gammas);
            else
                n = length(this.Distances);
            end
        end
        
        function applyConfiguration(this, nr, kernel)
            if ~isempty(this.Distances)
                this.Gammas(nr) = kernel.setGammaForDistance(this.Distances(nr),this.DistEps);
            else
                kernel.Gamma = this.Gammas(nr);
            end
        end
        
        function str = getConfigurationString(this, nr)
            str = [];
            if ~isempty(this.Gammas)
                str = sprintf('Gamma: %g',this.Gammas(nr));
                if ~isempty(this.Distances)
                    str = sprintf('%s (by dist %g)',str,this.Distances(nr));
                end
            end
        end
    end
    
end