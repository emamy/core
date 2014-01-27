classdef LogPlot
% LogPlot: Class with static functions for logarithmic plotting.
%
%
%
% @author Daniel Wirtz @date 2012-05-03
%
% @change{0,7,dw,2013-03-14} postprocessing does not call 'axis tight' if axis limits have
% already been set manually.
%
% @new{0,6,dw,2012-05-03} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods(Static)
        function p = nicesurf(h,X,Y,Z,varargin)
            % Creates a nice surface plot with the given data.
            %
            % Parameters:
            % h: The axes handle to plot to @type handle
            % X: The X data matrix @type matrix<double>
            % Y: The Y data matrix @type matrix<double>
            % Z: The `Z = \vf(X,Y)` value matrix @type matrix<double>
            % varargin: Any further arguments are passed to \c surf.
            %
            % Return values:
            % p: The patch object created by surf @type handle
            %
            % See also: surf
            p = surf(h,X,Y,Z,'FaceColor','interp','EdgeColor','k',varargin{:});
        end
        
        function p = nicesurfc(h,X,Y,Z,varargin)
            % Creates a nice surface plot including a contour with the given data.
            %
            % Parameters:
            % h: The axes handle to plot to @type handle
            % X: The X data matrix @type matrix<double>
            % Y: The Y data matrix @type matrix<double>
            % Z: The `Z = \vf(X,Y)` value matrix @type matrix<double>
            % varargin: Any further arguments are passed to \c surf.
            %
            % Return values:
            % p: The patch object created by surfc @type handle
            %
            % See also: surfc
            p = surfc(h,X,Y,Z,'FaceColor','interp','EdgeColor','k',varargin{:});
        end
        
        function p = logsurf(h,X,Y,Z,varargin)
            % Creates a nice surface plot in a logarithmic scale with the given data.
            %
            % Parameters:
            % h: The axes handle to plot to @type handle
            % X: The X data matrix @type matrix<double>
            % Y: The Y data matrix @type matrix<double>
            % Z: The `Z = \vf(X,Y)` value matrix @type matrix<double>
            % varargin: Any further arguments are passed to \c surf.
            %
            % Return values:
            % p: The patch object created by surf @type handle
            %
            % See also: surf
            
            % Create meshgrid if not already existing
            if isvector(X) && isvector(Y)
                [X,Y] = meshgrid(X,Y);
            end
            iszero = Z == 0;
            Z(iszero) = .5*min(Z(:));
            Z = log10(Z);
            p = LogPlot.nicesurf(h,X,Y,Z,varargin{:});
            LogPlot.postprocess(h);
        end
        
        function p = logsurfc(h,X,Y,Z,varargin)
            % Creates a nice surface plot including a contour with the given data in a
            % logarithmic scale.
            %
            % Parameters:
            % h: The axes handle to plot to @type handle
            % X: The X data matrix @type matrix<double>
            % Y: The Y data matrix @type matrix<double>
            % Z: The `Z = \vf(X,Y)` value matrix @type matrix<double>
            % varargin: Any further arguments are passed to \c surf.
            %
            % Return values:
            % p: The patch object created by surfc @type handle
            %
            % See also: surfc
            
            % Create meshgrid if not already existing
            if isvector(X) && isvector(Y)
                [X,Y] = meshgrid(X,Y);
            end
            iszero = Z == 0;
            Z(iszero) = .5*min(Z(:));
            Z = log10(Z);
            p = LogPlot.nicesurfc(h,X,Y,Z,varargin{:});
            LogPlot.postprocess(h);
        end
        
        function p = logtrisurf(h, tri, x, y, z, varargin)
            % Creates a surface plot from a 2D triangulation in a logarithmic scale.
            %
            % Parameters:
            % h: The axes handle to plot to @type handle
            % tri: The triangulation data, connecting the points given in the x,y variables.
            % x: The x values @type rowvec<double>
            % y: The y values @type rowvec<double>
            % z: The `z = \vf(x,y)` values @type rowvec<double>
            % varargin: Any further arguments are passed to \c trisurf.
            %
            % Return values:
            % p: The patch object(s) created by trisurf @type handle
            %
            % See also: trisurf
            z = log10(z);
            p = trisurf(tri, x, y, z,'Parent',h,'FaceColor','interp','EdgeColor','k',varargin{:});
            LogPlot.postprocess(h);
        end
        
        function p = cleverPlot(ax,x,y,varargin)
            % Calls corresponding plot routines depending on the scale of data.
            %
            % Parameters:
            % ax: The axes handle to plot to @type handle
            % x: The x values @type rowvec<double>
            % y: The y values @type rowvec<double>
            % varargin: Any further arguments are passed to the plot function.
            %
            % Return values:
            % p: The object(s) created by the plot function @type handle
            %
            % See also: plot semilogy
            if any(y(:)) < 0 || all(max(y)./min(y) < 50)
                pfun = @plot;
            else
                % As soon as one plot varies over more that two orders of magnitude,
                % use log Y scale
                set(ax,'YScale','log');
                pfun = @semilogy;
            end
            p = pfun(ax,x,y,varargin{:});
        end
    end
    
    methods(Static, Access=private)
        function postprocess(h)
            if all(strcmp('auto',{get(h,'XLimMode'),get(h,'YLimMode'),get(h,'ZLimMode')}))
                axis(h,'tight');
            end
            lbl = arrayfun(@(e)sprintf('%1.1e',e),10.^get(h,'ZTick'),'Unif',false);
            set(h,'ZTickLabel',lbl,'ZTickMode','manual');
        end
    end
end