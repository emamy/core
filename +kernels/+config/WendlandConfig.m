classdef WendlandConfig < kernels.config.RBFConfig
    % WendlandConfig: Configuration settings for Wendland kernels
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
        Smoothnesses;
    end
    
    properties
        Dimension;
    end
    
    methods
        function this = WendlandConfig(varargin)
            this = this@kernels.config.RBFConfig(varargin{:});
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('S',1);
            i.addParamValue('Dim',1,@(v)isscalar(v));
            i.parse(varargin{:});
            r = i.Results;
            this.Smoothnesses = r.S;
            if numel(this.Smoothnesses) ~= numel(this.Gammas)
                error('Need same amount of k values than gamma values');
            end
            this.Dimension = r.Dim;
            this.RequiredPrototypeClass = 'kernels.WendlandKernel';
        end
        
        function k = configureInstance(this, nr)
            k = configureInstance@kernels.config.RBFConfig(this, nr);
            k.k = this.Smoothnesses(nr);
            k.d = this.Dimension;
        end
        
        function str = getConfigurationString(this, nr, asCell)
            if nargin < 3
                asCell = false;
            end
            str = getConfigurationString@kernels.config.RBFConfig(this, nr, asCell);
            if asCell
                str = [str {sprintf('%d',this.Smoothnesses(nr))}]; %sprintf('%d',this.Dimension)
            else
                str = sprintf('%s, k: %d, d:%d',str,this.Smoothnesses(nr),this.Dimension);
            end
        end
        
        function str = getConfiguredPropertiesString(this)
            str = getConfiguredPropertiesString@kernels.config.RBFConfig(this);
            str = sprintf('%s, k with fixed d=%d',str,this.Dimension);
        end
        
        function conf = getSubPart(this, partNr, totalParts)
            conf = getSubPart@kernels.config.RBFConfig(this, partNr, totalParts);
            idx = this.getPartIndices(partNr, totalParts);
            conf.Smoothnesses = this.Smoothnesses(:,idx);
            conf.Dimension = this.Dimension;
        end
        
        function copy = clone(this)
            copy = kernels.config.WendlandConfig('G',[]);
            copy = clone@kernels.config.RBFConfig(this, copy);
            copy.Smoothnesses = this.Smoothnesses;
            copy.Dimension = this.Dimension;
        end
    end
    
    methods(Access=protected)
        function collectRanges(this, ptable, proppath)
            collectRanges@kernels.config.RBFConfig(this, ptable, proppath);
            this.addRange(ptable, [proppath {'Smoothnesses'}],min(this.Smoothnesses),max(this.Smoothnesses));
            this.addRange(ptable, [proppath {'Dimension'}],min(this.Dimension),max(this.Dimension));
        end
    end
end