classdef Polynomial < general.functions.AFunGen
    
    properties(Access=private)
        poly;
    end
    
    methods
        function this = Polynomial(poly)
            if nargin < 1
                poly = 1;
            end
            this.poly = poly;
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            p = this.poly;
            fhandle = @(t)polyval(p,t);
            dp = (8:-1:1) .* p(1:end-1);
            dfhandle = @(t)polyval(dp,t);
        end
        
        function str = getConfigStr(this)
            str = sprintf('Degree=%d',length(this.poly));
        end
    end
    
end