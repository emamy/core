classdef GaussConfig < kernels.config.RBFConfig
% RBFConfig: 
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
        Distances;
        DistEps;
    end
    
    methods
        function this = GaussConfig(varargin)
            this = this@kernels.config.RBFConfig(varargin{:});
            i = inputParser;
            i.KeepUnmatched = true;
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
                if isempty(this.Gammas)
                    error('You must pass either Gamma values (G) or distances (D)');
                end
            end
            this.RequiredPrototypeClass = 'kernels.GaussKernel';
        end
        
%         function str = getConfigurationString(this, nr, asCell)
%             if nargin < 3
%                 asCell = false;
%             end
%             str = getConfigurationString@kernels.config.RBFConfig(this, nr, asCell);
%             if ~asCell && ~isempty(this.Distances)
%                 str = sprintf('%s (by dist %g)',str,this.Distances(nr));
%             end
%         end
        
        function k = configureInstance(this, nr)
            if ~isempty(this.Distances)
                this.Gammas(nr) = kernel.setGammaForDistance(this.Distances(nr),this.DistEps);
            end
            k = configureInstance@kernels.config.RBFConfig(this, nr);
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = getSubPart@kernels.config.RBFConfig(this, partNr, totalParts);
            if ~isempty(this.Distances)
                idx = this.getPartIndices(partNr, totalParts);
                conf.Distances = this.Distances(:,idx);
                conf.DistEps = this.DistEps;
            end
        end
        
        function copy = clone(this)
            copy = kernels.config.GaussConfig('G',1);
            copy = clone@kernels.config.RBFConfig(this, copy);
            copy.Distances = this.Distances;
            copy.DistEps = this.DistEps;
        end
    end
end