classdef LogPlot
% LogPlot: Class with static functions for logarithmic plotting.
%
%
%
% @author Daniel Wirtz @date 2012-05-03
%
% @new{0,6,dw,2012-05-03} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        function p = logsurf(h,X,Y,Z,varargin)
            Z = log10(Z);
            p = surf(h,X,Y,Z,'FaceColor','interp','EdgeColor','interp',varargin{:});
            tools.LogPlot.postprocess(h);
        end
        
        function p = logtrisurf(h, tri, x, y, z, varargin)
            z = log10(z);
            p = trisurf(tri, x, y, z,'Parent',h,'FaceColor','interp','EdgeColor','interp',varargin{:});
            tools.LogPlot.postprocess(h);
            % color log-scaling mess
            %set(h,'zscale','log');
            %zdata = log10(get(p,'ZData'));
            %cdata = get(p,'CDATA');
        end
    end
    
    methods(Static, Access=private)
        function postprocess(h)
            axis(h,'tight');
            lbl = arrayfun(@(e)sprintf('%1.1e',e),10.^get(h,'ZTick'),'Unif',false);
            set(h,'ZTickLabel',lbl,'ZTickMode','manual');
        end
    end
end