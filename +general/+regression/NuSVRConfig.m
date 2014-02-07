classdef NuSVRConfig < IClassConfig
% NuSVRConfig: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-02-20
%
% @new{0,7,dw,2013-02-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(SetAccess=private)
        Nus = [];
        Lambdas = [];
    end
    
    methods
        function this = NuSVRConfig(values)
            if nargin == 1 && size(values,1) == 2
                this.Nus = values(1,:);
                this.Lambdas = values(2,:);
            end
            this.RequiredPrototypeClass = 'general.regression.ScalarNuSVR';
            this.Prototype = general.regression.ScalarNuSVR;
        end
        
        function n = getNumConfigurations(this)
            n = length(this.Nus);
        end
        
        function svr = configureInstance(this, nr)
            svr = this.getProtoClass;
            svr.nu = this.Nus(nr);
            svr.Lambda = this.Lambdas(nr);
        end
        
        function str = getConfigurationString(this, nr, ~)
            str = [];
            if ~isempty(this.Nus)
                str = sprintf('\nu: %g, \lambda: %g',this.Nus(nr),this.Lambdas(nr));
            end
        end
        
        function str = getConfiguredPropertiesString(~)
            str = 'nu, Lambda';
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = this.clone;
            idx = this.getPartIndices(partNr, totalParts);
            conf.Nus = this.Nus(idx);
            conf.Lambdas = this.Lambdas(idx);
        end
        
        function copy = clone(this)
            copy = general.regression.NuSVRConfig;
            copy = clone@IClassConfig(this, copy);
            copy.Nus = this.Nus;
            copy.Lambdas = this.Lambdas;
        end
        
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            this.addRange(ptable, [proppath {'nu'}],min(this.Nus),max(this.Nus));
            this.addRange(ptable, [proppath {'Lambda'}],min(this.Lambdas),max(this.Lambdas));
        end
    end
    
end