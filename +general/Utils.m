classdef Utils
% Collection of generally useful functions
%
% @author Daniel Wirtz @date 11.10.2010
%
% @change{0,6,dw,2012-06-08} Bugfix: Now the logarithmic norms are computed
% correctly. ARGH that took a while to detect..
%
% @new{0,6,dw,2012-04-13} New method "getTube" that allows to draw random
% vectors from spaces of arbitrary dimension but restricted to a specified
% tube.
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
            % glue: A string that is inserted between any element string
            % representation @type char @default ', '
            % format: The sprintf format string if data is a vector @type
            % char @default '%2.3e'
            %
            % Return values:
            % str: The concatented string of all 'data' strings glued together with the string
            % 'glue' @type char
            %
            % @new{0,6,dw,2011-11-17} Added the possibility to pass a numeric vector plus a
            % format string.
            str = '';
            if ~isempty(data)
                if nargin < 3
                    format = '%2.3e';
                    if nargin < 2
                        glue = ', ';
                    end
                end
                if isa(data,'cell')
                    str = data{1};
                    for idx = 2:length(data)
                        str = [str glue data{idx}];%#ok
                    end
                elseif isnumeric(data)
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
            
            if any(~ishandle(fig))
                error('Invalid figure handle.');
            end
            
            ExportDPI = '150';
            JPEGQuality = '95';
            exts = {'fig','pdf','eps','jpg','png','tif'};
            extd = {'MatLab Figure', 'PDF Files', 'Extended PostScript', 'JPEG Image',...
                'Portable Network Graphic', 'TIFF Image'};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            file = [];
            if nargin < 3
                extidx = 1;
                if nargin < 2
                    path = getpref('KERMOR','LASTPATH',pwd);
                    choices = cell(length(exts),2);
                    for i = 1:size(choices,1)
                        choices{i,1} = ['*.' exts{i}];
                        choices{i,2} = [extd{i} ' (*.' exts{i} ')'];
                    end
                    [filename, pathname, extidx] = uiputfile(choices, 'Save figure as', path);
                    if filename ~= 0 
                        file = [pathname filename];
                        setpref('KERMOR','LASTPATH',pathname);
                    end
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
%                 [d,fname] = fileparts(file);
                d = fileparts(file);
                if ~isempty(d) && exist(d,'file') ~= 7
                    mkdir(d);
                end
%                 allax = findobj(get(gcf,'Children'),'Type','axes');
%                 if length(allax) > 1
%                     if (extidx == 1)
%                         saveas(fig,file,'eps2c');
%                     elseif extidx == 2
%                         print(fig,file,['-djpeg' JPEGQuality],['-r' ExportDPI]);
%                     elseif extidx == 3
%                         saveas(fig, file, 'fig');
%                     elseif extidx == 4
%                         print(fig,file,'-dpdf',['-r' ExportDPI]);
%                     elseif extidx == 5
%                         saveas(fig,file,'png');
%                     end
%                 else
                    if extidx == 1 % fig
                        saveas(fig, file, 'fig');
                    else
                        args = {file, ['-' exts{extidx}],['-r' ExportDPI]};
                        if any(extidx == [2 3]) %pdf, eps
                            %args{end+1} = '-painters';
                            args{end+1} = '-transparent';
                        elseif extidx == 4 % jpg
                            args{end+1} = ['-q' JPEGQuality];
                            %args{end+1} = '-opengl';
                        elseif extidx == 5 % png
                            args{end+1} = '-transparent';
                        end
                        
%                         allax = findobj(get(fig,'Children'),'Type','axes');
%                         allax(strcmpi(get(allax,'Tag'),'legend')) = [];
%                         if length(allax) > 1
%                             args{end+1} = allax;
%                         else
                            args{end+1} = fig;
%                         end
                        
                        % export_fig ignores -transparent somehow..
                        c = get(fig,'Color');
                        set(fig,'Color','white');
                        
                        export_fig(args{:});
                        
                        set(fig,'Color',c);
                        %saveas(fig, fullfile(d,fname), 'fig');
                    end
%                 end
            else
                fprintf(2,'No file specified. Aborting\n');
            end
        end
        
        function saveAxes(ax, varargin)
            % Convenience function. Allows to save a custom axes instead of
            % a whole figure which allows to drop any unwanted uiobjects
            % contained on the source figure.
            %
            % @todo use export_fig here!
            
            fig = figure('Visible','off','MenuBar','none','ToolBar','none');
            %fig = figure('MenuBar','none','ToolBar','none');
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
            fu = get(f,'Units');
            au = get(a,'Units');
            
            set(f,'Units','pixels');
            set(a,'Units','pixels');
            fpos = get(f,'Position');
            apos = get(a,'Position');
            ati = get(a,'TightInset');
            set(f,'ActivePositionProperty','Position');
            set(a,'ActivePositionProperty','Position');
            set(f,'Position',[fpos(1:2) apos(3:4)+ati(1:2)+ati(3:4)]);
            set(a,'Position',[ati(1:2) apos(3:4)]);
            
            set(f,'Units',fu);
            set(a,'Units',au);
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
        
        function x = getTube(dim, num, length, spread, seed)
            % Computes `n=`'num' random vectors inside a tube of length
            % `l=`'length' inside a `d=`'dim' dimensional space.
            % The tube starts at zero and ends at a random point `x_e\in\R^d`
            % of length `l`.
            % The spread `s=`'spread' determines the maximum distance `ls`
            % by which the random vectors might be away from the tube
            % center `\{x\in\R^d~|~ x = \alpha x_e,\alpha\in[0,1]\}`.
            %
            % Parameters:
            % dim: The dimension of the tube. @type integer
            % num: The number of desired random vectors. @type integer
            % length: The length `l` of the tube. @type double @default 1
            % spread: The spread `s` of the tube vectors around the tube
            % core. @type double @default 0.15
            % seed: The seed to use for the random number generator. @type
            % integer @default cputime*10
            %
            % Return values:
            % x: An `d\times n` matrix with random vectors from a random
            % tube. @type matrix<double>
            if nargin < 5
                seed = cputime*10;
                if nargin < 4
                    spread = .15;
                    if nargin < 3
                        length = 1;
                    end
                end
            end
            r = RandStream('mt19937ar','Seed',seed);
            farend = r.rand(dim,1)-.5;
            farend = farend * length/norm(farend);
            x = r.rand(dim,num)-.5;
            x = bsxfun(@times, x, length*spread./Norm.L2(x));
            x = x + bsxfun(@times, farend, r.rand(1,num));
        end
        
        function [ln, v] = logNorm(A, G)
            % Computes the logarithmic norm of a matrix, optionally using a
            % positive definite symmetric matrix `G` inducing the matrix
            % norm to use.
            %
            % Parameters:
            % A: The target matrix @type matrix<double>
            % G: The positive definite, symmetric matrix `G` inducing the
            % norm and hence matrix norm to use. Defaults to identity
            % matrix (Euclidean/L2 norm) @type matrix<double>
            % @default I
            %
            % Return values:
            % ln: The logarithmic norm with respect to the given matrix
            % norm induced by `G`. @type double
            % v: The eigenvector of the largest eigenvalue (=log norm) of
            % the symmetric part of `A`. @type colvec<double>
            if nargin == 2
                L = chol(G,'lower');
                A = L\(A'*L);
            end

            hlp = .5*(A + A');
            if size(hlp,1) > 1000
%                 t = tic;
                opts.v0 = ones(size(A,1),1);
                opts.maxit = 1000;
                [v, ln] = eigs(hlp,1,'la',opts);
                if ln < 0
                    [v, ln] = eigs(hlp,1,-ln,opts);
%                     [v, lnpos] = eigs(hlp - ln*speye(size(A,1)),[],1,opts);
%                     ln = ln + lnpos;
                end
%                 t2 = toc(t);
            else
%                 t = tic;
                [v, d] = eig(full(hlp));
                [ln, idx] = max(diag(d));
                v = v(:,idx);
%                 t1 = toc(t);
            end
%             t = tic;
%             % Get smallest and largest magnitude eigenvalues
%             [v, ln2] = eigs(hlp,2,'be',opts);
%             [ln2, idx] = max([ln2(1) ln2(4)]);
%             v = v(:,idx);
%             t3 = toc(t);
%             fprintf('eig-full: %g/%gs, eigs-shift: %g/%gs, eigs-be: %g/%gs\n',...
%                 ln1,t1,ln,t2,ln2,t3);
        end
    end
    
    methods(Static)
        function res = test_Tube
            num = 1000;
            res = true;
            for i=1:50
                dim = round(rand*30);
                length = rand*100;
                spread = rand/2;
                
                x = general.Utils.getTube(dim, num, length, spread);
                res = res & all(Norm.L2(x) <= length*(1+spread));
            end
        end
        
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

