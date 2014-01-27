classdef ParamTimeExpansionConfig < kernels.config.ExpansionConfig
% ParamTimeExpansionConfig: Collects several class configurations for
% state, time and paramter kernels of a ParamTimeKernelExpansion
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
    
    properties
        % The time kernel configuration
        %
        % @type IClassConfig @default []
        TimeConfig = [];
        
        % The parameter kernel configuration
        %
        % @type IClassConfig @default []
        ParamConfig = [];
    end
    
    methods
        function this = ParamTimeExpansionConfig()
            this.RequiredPrototypeClass = 'kernels.ParamTimeKernelExpansion';
            this.Prototype = kernels.ParamTimeKernelExpansion;
        end
        
        function n = getNumConfigurations(this)
            % Returns the number of configurations that can be applied
            %
            % Return values:
            % n: The number of configurations @type integer
        
            % Passive strategy: Only use as many configurations as
            % available over all kernel configs.
            n = getNumConfigurations@kernels.config.ExpansionConfig(this);
            if ~isempty(this.TimeConfig)
                n = min(n,this.TimeConfig.getNumConfigurations);
            end
            if ~isempty(this.ParamConfig)
                n = min(n,this.ParamConfig.getNumConfigurations);
            end
        end
        
        function kexp = configureInstance(this, nr)
            kexp = configureInstance@kernels.config.ExpansionConfig(this, nr);
            if ~isempty(this.TimeConfig)
                kexp.TimeKernel = this.TimeConfig.configureInstance(nr);
            end
            if ~isempty(this.ParamConfig)
                kexp.ParamKernel = this.ParamConfig.configureInstance(nr);
            end
        end
        
        function str = getConfigurationString(this, nr, asCell)
            if nargin < 3
                asCell = false;
            end
            str = getConfigurationString@kernels.config.ExpansionConfig(this, nr, asCell);
            if ~isempty(this.TimeConfig)
                str = [str sprintf(', TimeKernel: %s',this.TimeConfig.getConfigurationString(nr, asCell))];
            end
            if ~isempty(this.ParamConfig)
                str = [str sprintf(', ParamKernel: %s',this.ParamConfig.getConfigurationString(nr, asCell))];
            end
        end
        
        function str = getConfiguredPropertiesString(this)
            e = {getConfiguredPropertiesString@kernels.config.ExpansionConfig(this)};
            if ~isempty(this.ParamConfig) || ~isempty(this.TimeConfig)
                e{end} = ['State: ' e{end}];
            end
            if ~isempty(this.TimeConfig)
                e(end+1) = {['Time: ' this.TimeConfig.getConfiguredPropertiesString]};
            end
            if ~isempty(this.ParamConfig)
                e(end+1) = {['Param: ' this.ParamConfig.getConfiguredPropertiesString]};
            end
            str = Utils.implode(e,', ');
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = kernels.config.ParamTimeExpansionConfig;
            if ~isempty(this.StateConfig)
                conf.StateConfig = this.StateConfig.getSubPart(partNr, totalParts);
            end
            if ~isempty(this.TimeConfig)
                conf.TimeConfig = this.TimeConfig.getSubPart(partNr, totalParts);
            end
            if ~isempty(this.ParamConfig)
                conf.ParamConfig = this.ParamConfig.getSubPart(partNr, totalParts);
            end
        end
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            collectRanges@kernels.config.ExpansionConfig(this, ptable, proppath);
            if ~isempty(this.TimeConfig)
                this.TimeConfig.collectRanges(ptable, ...
                    [proppath {sprintf('Time(%s)',this.TimeConfig.getClassName)}]);
            end
            if ~isempty(this.ParamConfig)
                this.ParamConfig.collectRanges(ptable, ...
                    [proppath {sprintf('Param(%s)',this.TimeConfig.getClassName)}]);
            end
        end
    end
    
    methods(Static)
        function res = test_ParamTimeExpansionConfig
            ptc = kernels.config.ParamTimeExpansionConfig;
            n = randi(100)*2;
            gk = kernels.config.GaussConfig('G',1:n);
            ptc.StateConfig = gk;
            ptc.ParamConfig = gk;
            tk = kernels.config.PolyConfig(linspace(0,5,n-1));
            ptc.TimeConfig = tk;
            
            res = ptc.getNumConfigurations == n-1;
            
            tk = kernels.config.PolyConfig(linspace(0,5,n));
            ptc.TimeConfig = tk;
            res = res && ptc.getNumConfigurations == n;
            
            kexp = ptc.configureInstance(1);
            res = res && isa(kexp,'kernels.ParamTimeKernelExpansion');
            
            disp(ptc.getConfigurationString(1));
            disp(ptc.getConfiguredPropertiesString);
            ptc1 = ptc.getSubPart(1,2);
            res = res && ptc1.getNumConfigurations == n/2;
            ptc2 = ptc.getSubPart(2,2);
            res = res && ptc2.getNumConfigurations == n/2;
            
            t = ptc.getValueRanges;
            t.print;
            
        end
    end
    
end