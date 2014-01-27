classdef ExpansionConfig < IClassConfig
% ExpansionConfig: Base class config for kernel expansions
%
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
        % The state space kernel configuration
        %
        % @type IClassConfig @default []
        StateConfig = [];
    end
    
    methods
        function this = ExpansionConfig()
            this.RequiredPrototypeClass = 'kernels.KernelExpansion';
            this.Prototype = kernels.KernelExpansion;
        end
        
        function n = getNumConfigurations(this)
            n = 0;
            if ~isempty(this.StateConfig)
                n = this.StateConfig.getNumConfigurations;
            end
        end
        
        function kexp = configureInstance(this, nr)
            kexp = this.getProtoClass;
            if ~isempty(this.StateConfig)
                kexp.Kernel = this.StateConfig.configureInstance(nr);
            end
        end
        
        function str = getConfigurationString(this, nr, asCell)
            if nargin < 3
                asCell = false;
            end
            str = [];
            if ~isempty(this.StateConfig)
                str = sprintf('StateKernel: %s',...
                    this.StateConfig.getConfigurationString(nr, asCell));
            end
        end
        
        function str = getConfiguredPropertiesString(this)
            str = this.StateConfig.getConfiguredPropertiesString;
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = kernels.config.ExpansionConfig;
            if ~isempty(this.StateConfig)
                conf.StateConfig = this.StateConfig.getSubPart(partNr, totalParts);
            end
        end
        
        function guessGammaConfigs(this, atd, ng, mf, Mf)
            % Guesses 'ng' `\gamma` values that might be geometrically senseful and sets the
            % applicable RBFConfig instances.
            %
            % Parameters:
            % atd: The approximation training data @type data.ApproxTrainData
            % ng: The number of `\gamma` values to guess @type integer
            % mf: The factor for minimum sample data distance @type double @default .5
            % Mf: The factor for the bounding box diameter @type double @default 2
            %
            
            if ~isa(atd,'data.ApproxTrainData')
                error('First argument has to be a data.ApproxTrainData instance');
            end
            gameps = .1;
            if nargin < 5
                Mf = 2;
                if nargin < 4
                    mf = .5;
                end
            end
            dfun = @logsp;
            
            dists = zeros(3,ng);
            dists(1,:) = dfun(mf*atd.xiDia, Mf*atd.xiDia);
            if atd.hasTime
                dists(2,:) = dfun(mf*atd.tiDia, Mf*atd.tiDia);
            end
            if atd.hasParams
                dists(3,:) = dfun(mf*atd.muiDia, Mf*atd.muiDia);
            end
            k = kernels.GaussKernel;
            this.Gammas = zeros(size(dists));
            for i=1:ng
               this.Gammas(1,i) = k.setGammaForDistance(dists(1,i),gameps);
               this.Gammas(2,i) = k.setGammaForDistance(dists(2,i),gameps);
               this.Gammas(3,i) = k.setGammaForDistance(dists(3,i),gameps);
            end
            
            function linsp(from, to)%#ok
                d = linspace(from,to,ng);
            end
            
            function d = logsp(from, to)
                d = logspace(log10(from),log10(to),ng);
            end
        end
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            if ~isempty(this.StateConfig)
                this.StateConfig.collectRanges(ptable, ...
                    [proppath {sprintf('State(%s)',this.StateConfig.getClassName)}]);
            end
        end
    end
    
end