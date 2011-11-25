classdef Utils
    % Collection of generally useful functions
    %
    % @author Daniel Wirtz @date 11.10.2010
    %
    % @change{0,6,dw,2011-11-17} Moved the getObjectConfig to a separate file.
    %
    % @change{0,6,dw,2011-11-16} Using mex CalcMD5 now to compute hash values for vectors.
    % Source downloaded from http://www.mathworks.com/matlabcentral/fileexchange/25921. Also
    % updated the KerMor.setup script to automatically compile the CalcMD5 mex file.
    %
    % @change{0,5,dw,2011-09-15} 
    % - saveAxes and saveFigure now store the last save location in the
    % preferences and reuse them.
    % - removeMargin now properly works, together with saveFigure or
    % saveAxes.
    %
    % @new{0,5,dw,2011-07-05} Added the general.Utils.implode function.
    %
    % @new{0,3,dw,2011-04-20} Added a new function general.Utils.getHelpShort to extract the first
    % line(s) of a help text in matlab style (text until first emtpy line = short)
    %
    % @new{0,3,dw,2011-04-18} Added the 'saveFigure' and 'saveAxes' methods from SegMedix.
    %
    % @change{0,3,dw,2011-04-04} Moved the general.Utils.getObjectConfig
    % method here from models.BaseModel
    %
    % @new{0,3,dw,2011-04-01}
    % - Added the general.Utils.getBoundingBox function.
    % - Added the general.Utils.findVecInMatrix function.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    methods(Static)
        
        function short = getHelpShort(txt)
            % Extracts the help short subtext from a given text.
            %
            % Gets the first block of a text that goes until the first
            % blank line.
            %
            % Parameters:
            % txt: The text to use. @type char
            %
            % Return values:
            % short: The subtext. @type char
            pos = regexp(txt,sprintf('\n[ ]*\n'));
            short = '';
            if ~isempty(pos)
                short = txt(1:pos(1)-1);
            else
                % Maybe only one line?
                pos = strfind(txt,char(10));
                if ~isempty(pos)
                    short = txt(1:pos(1)-1);
                end
            end
            short = strtrim(short);
        end
        
        function [bmin, bmax] = getBoundingBox(vectors)
            % Gets the bounding box for a matrix containing column vectors.
            % 
            % Parameters:
            % vectors: A `n\times m` matrix containing `m` column vectors
            % @type double
            % 
            % Return values:
            % bmin: A `n\times 1` vector representing the minimum value
            % corner of the bounding box @type double
            % bmax: A `n\times 1` vector representing the maximum value
            % corner of the bounding box @type double
            bmin = min(vectors,[],2);
            bmax = max(vectors,[],2);
        end
        
        function comb = createCombinations(ranges, varargin)
            % Creates the cartesian product of the vectors passed as a
            % matrix containing elements of each vector per row.
            %
            % Inputs:
            % ranges: Can either be a cell array of vectors or a vector.
            % varargin: If the first argument is a vector, an arbitrary
            % number of additional vectors can be passed to build the
            % cartesian product from.
            %
            % Return values:
            % comb: A matrix containing the combinations, each row
            % corresponds to an input vector's range.
            %
            % @author Daniel Wirtz @date 2010-10-11
            
            if ~isa(ranges,'cell')
                if isempty(varargin)
                    comb = ranges;
                    return;
                end
                r = cell(1,length(varargin)+1);
                r{1} = ranges;
                [r{2:end}] = varargin{:};
                ranges = r;
            end
            
            n = length(ranges);
            % Create nd-grids
            [matrices{1:n}] = ndgrid(ranges{:});
            % Convert to np x params matrix
            comb = zeros(n,numel(matrices{1}));
            for idx=1:n
                % Check if any range is empty - return empty then
                if isempty(matrices{idx})
                    comb = [];
                    return;
                end
                
                comb(idx,:) = matrices{idx}(:);
            end
        end
        
        function target = copyStructFields(source, target)
            % Recursively copies struct fields from one struct to another.
            %
            % Effectively implements a struct.clone() method.
            %
            % @author Daniel Wirtz @date 2010-11-03
            if ~isstruct(source) || ~isstruct(target)
                error('Both source and target arguments must be structs.');
            end
            % Get the field names from the source struct
            names = fieldnames(source);
            for idx = 1:length(names)
                % For struct fields, recursively copy the inner struct
                if isstruct(source.(names{idx}))
                    % Create target struct if not already set
                    if isempty(target.(names{idx}))
                        target.(names{idx}) = struct;
                    end
                    target.(names{idx}) = general.Utils.copyStructFields(source.(names{idx}),target.(names{idx}));
                % Else just copy the field values
                else
                    target.(names{idx}) = source.(names{idx});
                end
            end
        end
        
%         function str = getObjectConfig(obj, depth, numtabs)
%             % Gets a complete string representation of an object's
%             % configuration, sorted alphabetically.
%             %
%             % Parameters:
%             % obj: The object to get the configuration of
%             % depth: The maximum depth to go for
%             % sub-objects in properties. @type integer
%             % numtabs: [Optional, default 0] The number of tabs to insert
%             % before each output.
%             % Not necessary to set for normal calls as this is used upon
%             % recursive calls.
%             %
%             % Return values:
%             % str: The string representation of the object's state.
%             %
%             % @author Daniel Wirtz @date 2011-04-04
%             %
%             % @change{0,3,dw,2011-04-05} The order of the properties listed
%             % is now alphabetically, fixed no-tabs-bug.
%             %
%             % @todo remove from this class (currently sims are running who might call this)
%             str = object2str(obj);
%         end
        
        function str = implode(data, glue, format)
            % Implodes the elements of data using glue.
            %
            % Either transforms a cell array of strings into one string or implodes a numeric
            % vector using the specified format.
            %
            % Parameters:
            % data: A cell array of strings/chars or a numeric vector @type char|rowvec
            % glue: A string that is inserted between any element string representation @type char
            % format: The sprintf format string if data is a vector @type char
            %
            % Return values:
            % str: The concatented string of all 'data' strings glued together with the string
            % 'glue' @type char
            %
            % @new{0,6,dw,2011-11-17} Added the possibility to pass a numeric vector plus a
            % format string.
            str = '';
            if ~isempty(data)
                if isa(data,'cell')
                    str = data{1};
                    for idx = 2:length(data)
                        str = [str glue data{idx}];%#ok
                    end
                elseif isnumeric(data) && nargin == 3
                    % first n-1 entries
                    if numel(data) > 1
                        str = sprintf([format glue],data(1:end-1));
                    end
                    % append last, no glue afterwards needed
                    str = [str sprintf(format,data(end))];
                else
                    error('Can only pass cell arrays of strings or a vector with sprintf format pattern');
                end
            end
        end
        
        function idx = findVecInMatrix(A,b)
            % Finds column vectors inside a matrix.
            %
            % For multiple occurences, the first found index is used.
            %
            % See
            % http://www.mathworks.com/matlabcentral/newsreader/view_thread/174277
            % and the test function testing.find_vec_in_matrix_speedtest
            % for further information.
            %
            % Parameters:
            % A: A `n\times m` matrix of `m` column vectors
            % b: A `n\times p` vector, where each column is regarded as one vector to search
            %
            % Return values:
            % idx: A `1 \times p` vector containing the first found positions indices if a vector
            % from b is contained in A, zero otherwise.
            %
            % @change{0,3,dw,2011-04-12} Added support for multi vector search.
            % @change{0,3,dw,2011-04-13} Fixed errors when multiple occurences appear.
            
            if size(A,1) ~= size(b,1)
                error('Invalid arguments.');
            end
            idx = zeros(1,size(b,2));
            for n = 1:size(b,2)
                tmp = strfind(reshape(A,1,[]),b(:,n)');
                if ~isempty(tmp)
                    idx(n) = (tmp(1)+size(b,1)-1)/size(b,1);
                end
            end
        end
        
        function y = preparePlainPlot(y)
            % Memory-saving plotting for plain result plots.
            %
            % Parameters:
            % y: A result matrix `y` with rows corresponding to single
            % dimensions and columns corresponding to time-steps.
            %
            % Return values:
            % If there are more than 1000 dimensions, the subset with
            % distinct (via unique) last values are extracted. If this
            % still results in more than 4000 plots, the first 4000
            % dimensions are selected.
            if size(y,1) > 1000
                fprintf('Utils/preparePlainPlot: Number of graphs for plot > 1000, taking graphs with distinct y(T) values.\n');
                [v,idx] = unique(y(:,end));
                [v,idxm] = unique(y(:,round(size(y,2)/2)));
                y = y(union(idx,idxm),:);
                sy = size(y,1);
                if sy > 4000
                    fprintf('Utils/preparePlainPlot: Number of graphs for plot still > 4000, taking 4000 graphs.\n');
                    y = y(round(1:sy/4000:sy),:);
                end
            end
        end
        
        function saveFigure(fig, filename, ext)
            % Opens a matlab save dialog and saves the given figure to the
            % file selected.
            %
            % Supported formats: eps, jpg, fig
            %
            % @change{0,4,dw,2011-05-31} Improved the export capabilites and automatic removement of
            % any figure margins is performed.
            
            ExportDPI = '300';
            JPEGQuality = '95';
            exts = {'eps','jpg','fig','pdf','png'};
            
            if nargin < 3
                extidx = 1;
                if nargin < 2
                    path = getpref('KERMOR','LASTPATH',pwd);
                    [filename, pathname, extidx] = ...
                        uiputfile({'*.eps','Extended PostScript (*.eps)';...
                       '*.jpg','JPEG Image (*.jpg)';...
                       '*.fig','Figures (*.fig)';...
                       '*.pdf','PDF Files (*.pdf)';...
                       '*.png','Portable Network Graphic (*.png)}'}, 'Save figure as',path);
                    file = [pathname filename];
                    setpref('KERMOR','LASTPATH',pathname);
                else
                    file = [filename '.' exts{extidx}];
                end
            else
                extidx = find(strcmp(ext,exts),1);
                if isempty(extidx)
                    warning('KerMor:Utils:invalidExtension','Invalid extension: %s, using eps',ext);
                    extidx = 1;
                end
                file = [filename '.' exts{extidx}];
            end
            
            if ~isempty(file)
                if length(get(gcf,'Children')) > 1
                    if (extidx == 1)
                        saveas(fig,file,'eps2c');
                    elseif extidx == 2
                        print(fig,file,['-djpeg' JPEGQuality],['-r' ExportDPI]);
                    elseif extidx == 3
                        saveas(fig, file, 'fig');
                    elseif extidx == 4
                        print(fig,file,'-dpdf',['-r' ExportDPI]);
                    elseif extidx == 5
                        saveas(fig,file,'png');
                    end
                else
                    a = gca(fig);

                    % Store old position and remove margin for export
                    oldap = get(a,'ActivePosition');
                    oldpos = get(a,'Position');
                    %oldmode = get(f,'PaperPositionMode');
                    set(fig,'PaperPositionMode','auto');%,'InvertHardcopy','off');

                    if (extidx == 1)
                        %print(fig,file,'-depsc2',['-r' ExportDPI],'-tiff');
                        general.Utils.removeMargin(fig);
                        saveas(fig,file,'eps2c');
                        %system(['xdg-open ' file]);
                    elseif extidx == 2
                        general.Utils.removeMargin(fig);
                        print(fig,file,['-djpeg' JPEGQuality],['-r' ExportDPI]);
                        %system(['xdg-open ' file]);
                    elseif extidx == 3
                        saveas(fig, file, 'fig');
                        %openfig(file,'new','visible');
                    elseif extidx == 4
                        general.Utils.removeMargin(fig);
                        print(fig,file,'-dpdf',['-r' ExportDPI]);
                    elseif extidx == 5
                        general.Utils.removeMargin(fig);
                        saveas(fig,file,'png');
                    end

                    % Restore old position
                    set(a,'ActivePosition',oldap);
                    set(a,'Position',oldpos);
                end
            end
        end
        
        function saveAxes(ax, varargin)
            % Convenience function. Allows to save a custom axes instead of
            % a whole figure which allows to drop any unwanted uiobjects
            % contained on the source figure.
            
            fig = figure('Visible','off','MenuBar','none','ToolBar','none');
            %fig = figure;
            % Set fig size to axis size
            %set(fig,'Position', [fpos(1:2) (apos(3:4)+ti(3:4))/2]);
            copyobj(ax, fig);
            
            %% Fit style
            % Just copy the colormap
            %set(fig,'Colormap',get(get(ax,'Parent'),'Colormap'));
                        
            %% Save
            general.Utils.saveFigure(fig, varargin{:});
            close(fig);
        end
        
        function removeMargin(f)
            % Requires the axes and figure units to be the same.
            %
            % @todo tailor this method so that subplots are also supported (so far only one
            % axis!) and change saveFigure again.
            a = gca(f);
            set(f,'Units','pixels');
            set(a,'Units','pixels');
            fpos = get(f,'Position');
            apos = get(a,'Position');
            ati = get(a,'TightInset');
            set(f,'ActivePositionProperty','Position');
            set(a,'ActivePositionProperty','Position');
            set(f,'Position',[fpos(1:2) apos(3:4)+ati(1:2)+ati(3:4)]);
            set(a,'Position',[ati(1:2) apos(3:4)]);
        end
        
        function h = getHash(vec)
            % Returns a hash code for the given vector.
            %
            % Currently using the CalcMD5 routine from 3rdparty/calcmd5, which is included in
            % KerMor as 3rd party code.
            % Original source: http://www.mathworks.com/matlabcentral/fileexchange/25921
            % 
            h = CalcMD5(vec);
            %h = sprintf('%d',typecast(vec,'uint8'));
        end
        
    end
    
    methods(Static)
        function res = test_createCombinations
            % Tests the createCombinations function.
            % @author Daniel Wirtz @date 11.10.2010
            res = true;
            
            res = res && isequal([1 2 3 1 2 3; 1 1 1 2 2 2],general.Utils.createCombinations(1:3,1:2));
            
            res = res && isempty(general.Utils.createCombinations(1:3,1:2,[],1:54));
            
            res = res && isequal(1:20,general.Utils.createCombinations(1:20));
        end
        
        function res = test_findVec
            % Tests the findVecInMatrix function.
            % @author Daniel Wirtz @date 2011-04-12
            a = rand(40,40);
            idx = randperm(40);
            idx = idx(1:10);
            b = a(:,idx);
            
            % Add nonexistent
            b(:,end+1) = rand(40,1);
            idx = [idx 0];
            
            % Add multiples
            a(:,end+1:end+2) = a(:,[idx(1) idx(5)]);
            
            res = all(idx - general.Utils.findVecInMatrix(a,b) == 0);
        end
    end
    
end

