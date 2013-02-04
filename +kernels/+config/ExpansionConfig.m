classdef ExpansionConfig < general.IClassConfig
% ExpansionConfig: Collects several IKernelConfigs to apply to a kernel expansion
%
% @todo implement setter etc and do checks
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
        StateConfig = [];
        TimeConfig = [];
        ParamConfig = [];
    end
    
    methods
        function n = getNumConfigurations(this)
            n = Inf;
            if ~isempty(this.StateConfig)
                n = min(n,this.StateConfig.getNumConfigurations);
            end
            if ~isempty(this.TimeConfig)
                n = min(n,this.TimeConfig.getNumConfigurations);
            end
            if ~isempty(this.ParamConfig)
                n = min(n,this.ParamConfig.getNumConfigurations);
            end
            if isinf(n)
                n = 0;
            end
        end
        
        function applyConfiguration(this, nr, expansion)
            % @todo improve code (dont think str as return value for maybe verbose is any good)
            if ~isempty(this.StateConfig)
                this.StateConfig.applyConfiguration(nr, expansion.Kernel);
            end
            if ~isempty(this.TimeConfig)
                this.TimeConfig.applyConfiguration(nr, expansion.TimeKernel);
            end
            if ~isempty(this.ParamConfig)
                this.ParamConfig.applyConfiguration(nr, expansion.ParamKernel);
            end
        end
        
        function str = getConfigurationString(this, nr, asCell)
            if nargin < 3
                asCell = false;
            end
            str = [];
            if ~isempty(this.StateConfig)
                str = [str sprintf('state: %s',this.StateConfig.getConfigurationString(nr, asCell))];
            end
            if ~isempty(this.TimeConfig)
                str = [str sprintf('time: %s',this.TimeConfig.getConfigurationString(nr, asCell))];
            end
            if ~isempty(this.ParamConfig)
                str = [str sprintf('param: %s',this.ParamConfig.getConfigurationString(nr, asCell))];
            end
        end
        
        function str = getConfiguredPropertiesString(this)
            e = {};
            if ~isempty(this.StateConfig)
                if isempty(this.ParamConfig) && isempty(this.TimeConfig)
                    e(end+1) = {this.StateConfig.getConfiguredPropertiesString};
                else
                    e(end+1) = {['State: ' this.StateConfig.getConfiguredPropertiesString]};
                end
            end
            if ~isempty(this.TimeConfig)
                e(end+1) = {['Time: ' this.TimeConfig.getConfiguredPropertiesString]};
            end
            if ~isempty(this.ParamConfig)
                e(end+1) = {['Param: ' this.ParamConfig.getConfiguredPropertiesString]};
            end
            str = general.Utils.implode(e,', ');
        end
        
        function setBestConfig(this, idx, expansion)
            % @todo implement event listener for PostSet of vBestConfigIndex
            
            % Apply config
            this.applyConfiguration(idx, expansion);
            % Set best config indices
            this.vBestConfigIndex = idx;
            this.StateConfig.vBestConfigIndex = idx;
            if ~isempty(this.TimeConfig)
                this.TimeConfig.vBestConfigIndex = idx;
            end
            if ~isempty(this.ParamConfig)
                this.ParamConfig.vBestConfigIndex = idx;
            end
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = kernels.config.ExpansionConfig;
            if ~isempty(this.StateConfig)
                conf.StateConfig = this.StateConfig.getSubPart(this, partNr, totalParts);
            end
            if ~isempty(this.TimeConfig)
                conf.TimeConfig = this.TimeConfig.getSubPart(this, partNr, totalParts);
            end
            if ~isempty(this.ParamConfig)
                conf.ParamConfig = this.ParamConfig.getSubPart(this, partNr, totalParts);
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
    
end