classdef LineSpecIterator < handle
% LinSpecIterator: Small helper class to automatically iterate through different line
% styles/markers/colors
%
% @author Daniel Wirtz @date 2012-06-08
%
% @new{0,6,dw,2012-07-26} Can now also iterate through colors with the option of excluding
% some.
%
% @new{0,6,dw,2012-06-08} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties
        MarkerStyles = {'s','d','o','*','p','x','+','^','v','<','>','.','h'};
        
        LineStyles = {'-',':','-.','--'};
        
        Colors = [  0       0       1.0000;
                    0       0.5000  0;
                    1.0000  0       0;
                    0       0.7500  0.7500;
                    0.7500  0       0.7500;
                    0.7500  0.7500  0;
                    0.2500  0.2500  0.2500];
    end
    
    properties(Access=private)
        curl = 0;
        curm = 0;
        curc = 0;
        excluded_cols = [];
    end
    
    methods
        function markerstyle = nextMarkerStyle(this)
            markerstyle = this.MarkerStyles{this.curm+1};
            this.curm = mod(this.curm+1,length(this.MarkerStyles));
        end
        
        function linestyle = nextLineStyle(this)
            linestyle = this.LineStyles{this.curl+1};
            this.curl = mod(this.curl+1,length(this.LineStyles));
        end
        
        function color = nextColor(this)
            color = this.Colors(this.curc+1,:);
            this.curc = mod(this.curc+1,size(this.Colors,1));
            cnt = 1;
            while any(sum(abs(repmat(color,size(this.excluded_cols,1),1)-this.excluded_cols),2) == 0)
                % Pick next one
                color = this.Colors(this.curc+1,:);
                this.curc = mod(this.curc+1,size(this.Colors,1));
                
                cnt = cnt+1;
                if cnt == size(this.Colors,1)
                    warning('KerMor:LineSpecInterator',...
                        'Out of pre-defined colors (all excluded). Using random color.');
                    color = rand(1,3);
                    this.Colors = [this.Colors; color];
                    break;
                end
            end
        end
        
        function excludeColor(this, data)
            if size(data,2) == 3
                color = data;
            elseif ishandle(data)
                color = [];
                for i=1:length(data)
                    color = [color; get(data(i),'Color')];%#ok
                end
            end
            this.excluded_cols = [this.excluded_cols; color];
        end
    end
    
end