classdef EpsSVRConfig < IClassConfig
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetAccess=private)
        Epsilons = [];
        Lambdas = [];
    end
    
    methods
        function this = EpsSVRConfig(values)
            if nargin == 1 && size(values,1) == 2
                this.Epsilons = values(1,:);
                this.Lambdas = values(2,:);
            end
            this.RequiredPrototypeClass = 'general.regression.BaseScalarSVR';
            this.Prototype = general.regression.ScalarEpsSVR_SMO;
        end
        
        function n = getNumConfigurations(this)
            n = length(this.Epsilons);
        end
        
        function svr = configureInstance(this, nr)
            svr = this.getProtoClass;
            svr.Eps = this.Epsilons(nr);
            svr.Lambda = this.Lambdas(nr);
        end
        
        function str = getConfigurationString(this, nr, ~)
            str = [];
            if ~isempty(this.Epsilons)
                str = sprintf('\\epsilon: %g, \\lambda: %g',this.Epsilons(nr),this.Lambdas(nr));
            end
        end
        
        function str = getConfiguredPropertiesString(~)
            str = 'Eps, Lambda';
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = this.clone;
            idx = this.getPartIndices(partNr, totalParts);
            conf.Epsilons = this.Epsilons(idx);
            conf.Lambdas = this.Lambdas(idx);
        end
        
        function copy = clone(this)
            copy = general.regression.EpsSVRConfig;
            copy = clone@IClassConfig(this, copy);
            copy.Epsilons = this.Epsilons;
            copy.Lambdas = this.Lambdas;
        end
        
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            this.addRange(ptable, [proppath {'Eps'}],min(this.Epsilons),max(this.Epsilons));
            this.addRange(ptable, [proppath {'Lambda'}],min(this.Lambdas),max(this.Lambdas));
        end
    end
    
end