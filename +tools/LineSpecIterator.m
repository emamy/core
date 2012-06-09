classdef LineSpecIterator < handle
% LinSpecIterator: 
%
%
%
% @author Daniel Wirtz @date 2012-06-08
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
    end
    
    properties(Access=private)
        curl = 0;
        curm = 0;
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
    end
    
end