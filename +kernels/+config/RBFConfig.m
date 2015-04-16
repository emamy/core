classdef RBFConfig < IClassConfig
% RBFConfig: Base configuration settings for kernels implementing ARBFKernel.
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
    
    properties(SetAccess=protected)
        Gammas;
    end
    
    methods
        function this = RBFConfig(varargin)
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('G',[]);
            i.parse(varargin{:});
            r = i.Results;
            if isempty(r.G)
                this.Gammas = 1;
            else
                this.Gammas = r.G;
            end
            this.RequiredPrototypeClass = 'kernels.ARBFKernel';
        end
        
        function n = getNumConfigurations(this)
            n = length(this.Gammas);
        end
        
        function k = configureInstance(this, nr)
            k = this.getProtoClass;
            k.Gamma = this.Gammas(nr);
        end

        function str = getConfigurationString(this, nr, asCell)
            if nargin < 3
                asCell = false;
            end
            if asCell
                str = {};
                if ~isempty(this.Gammas)
                    str{1} = sprintf('%g',this.Gammas(nr));
                end
            else
                str = [];
                if ~isempty(this.Gammas)
                    str = sprintf('\\gamma: %g',this.Gammas(nr));
                end
            end
        end
        
        function str = getConfiguredPropertiesString(~)
            str = '\gamma';
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = this.clone;
            conf.Gammas = this.Gammas(this.getPartIndices(partNr, totalParts));
        end
        
        function copy = clone(this, copy)
            if nargin < 2
                copy = kernels.config.RBFConfig;
            end
            copy = clone@IClassConfig(this, copy);
            copy.Gammas = this.Gammas;
        end
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            this.addRange(ptable, [proppath {'gamma'}],min(this.Gammas),max(this.Gammas));
        end
    end
    
    methods(Static)
        function dists = getDists(atd, num)
            % Computes the distances for the different RBF kernel
            % Gamma configurations using the 'atd' data and the algorithms
            % configuration.
            %
            % Parameters:
            % atd: The approximation training data @type data.ApproxTrainData
            % num: The number of gamma values to compute @type integer
            %
            % Return values:
            % dists: A `3\times n` matrix with distance values for state,
            % time and parameter kernels (time and parameter if given, but
            % always 2nd and 3rd rows, respectively)
            
            dfun = @logsp; % distances comp fun (linsp / logsp)
            
            Mf = 2;
            mf = .5;
            if isscalar(mf)
                mf = ones(3,1)*mf;
            end
            if isscalar(Mf)
                Mf = ones(3,1)*Mf;
            end
            xm = atd.xiDia;
            dists = dfun(mf(1)*xm, Mf(1)*atd.xiDia);
            if atd.hasTime
                tm = atd.tiDia;
                dists(2,:) = dfun(mf(2)*tm, Mf(2)*atd.tiDia);
            end
            if atd.hasParams
                mum = atd.muiDia;
                dists(3,:) = dfun(mf(3)*mum, Mf(3)*atd.muiDia);
            end
            
            function d = linsp(from, to)%#ok
                d = linspace(from,to,num);
            end
            
            function d = logsp(from, to)
                d = logspace(log10(from),log10(to),num);
            end
        end
    end
    
end