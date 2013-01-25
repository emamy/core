classdef RBFConfig < general.IClassConfig
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
            this.Gammas = r.G;
        end
        
        function n = getNumConfigurations(this)
            if ~isempty(this.Gammas)
                n = length(this.Gammas);
            else
                n = length(this.Distances);
            end
        end
        
        function applyConfiguration(this, nr, kernel)
            kernel.Gamma = this.Gammas(nr);
        end
        
        function str = getConfigurationString(this, nr, asCell)
            if asCell
                str = {};
                if ~isempty(this.Gammas)
                    str{1} = sprintf('%g',this.Gammas(nr));
                end
            else
                str = [];
                if ~isempty(this.Gammas)
                    str = sprintf('Gamma: %g',this.Gammas(nr));
                end
            end
        end
        
        function str = getConfiguredPropertiesString(~)
            str = 'Gamma';
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