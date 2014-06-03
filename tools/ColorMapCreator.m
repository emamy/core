classdef ColorMapCreator < handle
% ColorMapCreator: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-06-02
%
% @new{0,7,dw,2014-06-02} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(Constant)
        Resolution = 1000;
    end
    
    properties
        LogPlot = false;
    end
    
    properties(Access=private)
        vals = [];
        cols = {};
    end
    
    methods
        
        function useJet(this, vals)
            this.addColor(vals(1),[0 0 1]);
            this.addColor(vals(2),[0 1 1]);
            this.addColor(vals(3),[1 1 0]);
            this.addColor(vals(4),[1 0 0]);
        end
        
        function addColor(this, value, color, colorbefore)
            if nargin < 4
                colorbefore = color;
            end
            this.vals(end+1) = value;
            this.cols{end+1} = color;
            if ~isequal(color, colorbefore)
                this.vals(end+1) = value-1/this.Resolution;
                this.cols{end+1} = colorbefore;
            end
        end
        
        function cmap = create(this, Zdata)
            dvals = this.vals;
            dcols = this.cols;
            
            if this.LogPlot
                dvals = log10(dvals);
            end
            
            if nargin > 1
                if this.LogPlot
                    Zdata = log10(Zdata);
                end
                minZ = min(Zdata(:));
                dvals(end+1) = minZ;
                dcols{end+1} = [0 0 .5]; %zeros(3,1)
                
                maxZ = max(Zdata(:));
                dvals(end+1) = maxZ;
                dcols{end+1} = [.5 0 0]; %zeros(3,1);
                
                rem = dvals < minZ;
                dvals(rem)  = [];
                dcols(rem) = [];
                rem = dvals > maxZ;
                dvals(rem)  = [];
                dcols(rem) = [];
            end
            
            % Colormap
            cmap = [];
            if length(dvals) > 1
                [val, idx] = sort(dvals,'ascend');
                valrange = linspace(val(1),val(end),this.Resolution);
                ncol = length(valrange);
                r = zeros(ncol,1);
                g = zeros(ncol,1);
                b = zeros(ncol,1);

                curpos = 1;
                curcol = dcols{idx(1)};
                for k = 2:length(val)
                    nextcol = dcols{idx(k)};
                    [~,nextpos] = min(abs(valrange-val(k)));

                    pos = curpos:nextpos;
                    r(pos) = linspace(curcol(1),nextcol(1),length(pos));
                    g(pos) = linspace(curcol(2),nextcol(2),length(pos));
                    b(pos) = linspace(curcol(3),nextcol(3),length(pos));

                    curpos = nextpos;
                    curcol = nextcol;
                end
                cmap = [r g b];
            end
        end
    end
    
end