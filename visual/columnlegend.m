function [lh, object_h, varargout] = columnlegend(varargin)
% legend: Overrides the legend builtin function and adds the index of the associated plot to
% the created line/text objects.
%
% @author Daniel Wirtz @date 2012-07-26
%
% @new{0,6,dw,2012-07-26} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

error('Not yet finished');

p = inputParser;
p.CaseSensitive = false;
p.addOptional('EnsureColumns',[],@(l)ishandle(l) && strcmp(get(l,'type'),'axes') && strcmp(get(l,'tag'),'legend'));
p.addOptional('Columns', 1, @(c)~isempty(c) && isscalar(c) && c>0 && round(c) == c);
p.addOptional('MaxRows', 5, @(r)~isempty(r) && isscalar(r) && r>0 && round(r) == r);
p.parse(varargin{:});
if ~isempty(p.Results.EnsureColumns)
    lh = p.Results.EnsureColumns;
    ud = get(lh,'UserData');
    
    columns = p.Results.Columns;
    
    numentries = length(get(l,'String'));
    
    applyColumnLayout(lh);
    
else % "normal" legend creation + possible layout
    [lh, object_h, varargout{1:nargout-2}] = legend(varargin{:});

    numentries = length(object_h)/3;
    for entrynr=1:numentries
        % Text
        text = object_h(entrynr);
        if ~strcmpi(get(text,'type'),'text')
            warning('KerMor:indexedlegend','Element was expected to be of type "text", found %s instead.',get(text,'type'));
        end
        set(text,'UserData',entrynr);

        % Line
        line = object_h(numentries + 2*(entrynr-1)+1);
        if ~strcmpi(get(line,'type'),'line')
            warning('KerMor:indexedlegend','Element was expected to be of type "line", found %s instead.',get(line,'type'));
        elseif size(get(line,'xdata'),2) ~= 2
            warning('KerMor:indexedlegend','XData does not indicate this element to be a non-point line object');
        end
        set(line,'UserData',entrynr);

        % Marker
        marker = object_h(numentries + 2*(entrynr-1)+2);
        if ~strcmpi(get(line,'type'),'line')
            warning('KerMor:indexedlegend','Element was expected to be of type "line", found %s instead.',get(marker,'type'));
        elseif size(get(marker,'xdata'),2) ~= 1
            warning('KerMor:indexedlegend','XData does not indicate this element to be marker (one point line) object');
        end
        set(marker,'UserData',entrynr);
    end
    
    % Only rearrange if it makes sense
    if (numentries > 1 && p.Results.Columns > 1) || numentries > p.Results.MaxRows
        
        % Default
        columns = p.Results.Columns;
        
        
        columns = ceil(numentries/this.MaxLegendRows);
    end
end

    function applyColumnLayout(lh, numentries, cols)

        % Only fuzz around if more than one column
        if columns > 1

            location = get(leg,'Location');
            % Set default location
            if strcmpi(location,'none')
                location = 'SouthEast';
                set(leg,'Location',location);
            end

            % get old width, new width and scale factor
            pos = get(leg, 'position');
            width = columns*pos(3);
            rescale = pos(3)/width;

            childs = get(leg,'Children');
            % get some old values so we can scale everything later
            xdata = get(childs(2), 'xdata'); % any x line position will do
            % Pick the this.MaxLegendRows' element from top (corresponds to the
            % MaxLegendRows' element from bottom in legend) to get the correct new y
            % height!
            ydata1 = get(childs(3*(this.MaxLegendRows-1)+1), 'ydata');
            % Take the element beneath that (surely exists)
            ydata2 = get(childs(3*(this.MaxLegendRows-1)-2), 'ydata');

            %we'll use these later to align things appropriately
            sheight = abs(ydata1(1)-ydata2(1));                  % height between data lines
            height = ydata1(1);                             % height of the box. Used to top margin offset
            line_width = (xdata(2)-xdata(1))*rescale;   % rescaled linewidth to match original
            spacer = xdata(1)*rescale;                    % rescaled spacer used for margins

            % put the legend on the lower left corner to make initial adjustments easier
            loci = get(gca, 'position');
            set(leg, 'position', [loci(1) pos(2) width pos(4)]);

            col = -1;
            for i=1:rows
                if mod(i,this.MaxLegendRows)==1,
                    col = col+1;
                end

                position = mod(i,this.MaxLegendRows);
                if position == 0,
                     position = this.MaxLegendRows;
                end

                elem = 3*rows-2 - (i-1)*3; % element index (last to first)
                ypos = height-(position-1)*sheight;
                % Marker
                set(childs(elem), 'xdata', col/columns+spacer*3.5,'ydata', ypos);

                % Lines
                set(childs(elem+1), ...
                    'xdata', [col/columns+spacer col/columns+spacer+line_width],...
                    'ydata', [ypos ypos]);

                % Labels
                set(childs(elem+2), 'position', [col/columns+spacer*2+line_width ypos]);
            end

            %unfortunately, it is not possible to force the box to be smaller than the
            %original height, therefore, turn it off and set background color to none
            %so that it no longer appears
            set(leg, 'Color', 'None', 'Box', 'off');

            % Get updated position
            pos = get(leg, 'position');
            % Get figure position
            ax_pos = get(gca, 'position');
            switch lower(location),
                case {'northeast'}
                    set(leg, 'position', [pos(1)+ax_pos(3)-pos(3) pos(2) pos(3) pos(4)]);
                case {'southeast'}
                    set(leg, 'position', [pos(1)+ax_pos(3)-pos(3) ax_pos(2) pos(3) pos(4)]);
                case {'southwest'}
                    set(leg, 'position', [ax_pos(1) ax_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
                otherwise
                    warning('KerMor:PlotManager','Cannot handle location ''%s'' correctly.',location);
            end
        end
    end
    

end