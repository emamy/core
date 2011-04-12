classdef Utils
    % Collection of generally useful functions
    %
    % @author Daniel Wirtz @date 11.10.2010
    %
    % @change{0,3,dw,2011-04-04} Moved the general.Utils.getObjectConfig
    % method here from models.BaseModel
    %
    % @new{0,3,dw,2011-04-01}
    % - Added the general.Utils.getBoundingBox function.
    % - Added the general.Utils.findVecInMatrix function.
    
    methods(Static)
        
        function idx = findVecInMatrix(A,b)
            % Finds a column vector inside a matrix.
            %
            % See
            % http://www.mathworks.com/matlabcentral/newsreader/view_thread/174277
            % and the test function testing.find_vec_in_matrix_speedtest
            % for further information.
            %
            % Parameters:
            % A: A `n\times m` matrix of `m` column vectors
            % b: A `n\times 1` vector
            %
            % Return values:
            % idx: The column index idx `\in \{1 \ldots m\}` if the vector
            % b is contained in A, [] otherwise.
            if size(b,2) ~= 1 || size(A,1) ~= size(b,1)
                error('Invalid arguments.');
            end
            idx = strfind(reshape(A,1,[]),b');
            idx = (idx+length(b)-1)/length(b);
        end
        
        function [bmin, bmax] = getBoundingBox(vectors)
            % Gets the bounding box for a matrix containing column vectors.
            % 
            % Parameters:
            % vectors: A `n\times m` matrix containing `m` column vectors
            % 
            % Return values:
            % bmin: A `n\times 1` vector representing the minimum value
            % corner of the bounding box
            % bmax: A `n\times 1` vector representing the maximum value
            % corner of the bounding box
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
            % @author Daniel Wirtz @date 11.10.2010
            
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
            % @author Daniel Wirtz @date 03.11.2010
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
        
        function str = getObjectConfig(obj, depth, numtabs)
            % Gets a complete string representation of an object's
            % configuration, sorted alphabetically.
            %
            % Parameters:
            % obj: The object to get the configuration of
            % depth: [Optional, default 5] The maximum depth to go for
            % sub-objects in properties.
            % numtabs: [Optional, default 0] The number of tabs to insert
            % before each output.
            % Not necessary to set for normal calls as this is used upon
            % recursive calls.
            %
            % Return values:
            % str: The string representation of the object's state.
            %
            % @author Daniel Wirtz @date 2011-04-04
            %
            % @change{0,3,dw,2011-04-05} The order of the properties listed
            % is now alphabetically, fixed no-tabs-bug.
            str = '';
            if nargin < 3
                numtabs = 0;
                if nargin < 2
                    depth = 5;
                elseif depth == 0
                    warning('Un:important','Exceeded the maximal recursion depth. Ignoring deeper objects.');
                    return;
                end
            end
            mc = metaclass(obj);
            if isfield(obj,'Name')
                name = obj.Name;
            else
                name = mc.Name;
            end
            if ~isempty(str)
                str = [str ': ' name '\n'];
            else
                str = [name ':\n'];
            end
            % get string cell of names and sort alphabetically
            names = cellfun(@(mp)mp.Name,mc.Properties,'UniformOutput',false);
            [val, sortedidx] = sort(names);
            for n = 1:length(sortedidx)
                idx = sortedidx(n);
                p = mc.Properties{idx};
                if strcmp(p.GetAccess,'public')
                    str = [str repmat('\t',1,numtabs)];%#ok
                    pobj = obj.(p.Name);
                    if ~isempty(pobj) && ~isequal(obj,pobj)
                        if isobject(pobj)
                            str = [str p.Name general.Utils.getObjectConfig(pobj, depth-1, numtabs+1)];%#ok
                        elseif isnumeric(pobj)
                            if numel(pobj) > 20
                                str = [str p.Name ': [' num2str(size(pobj)) '] ' class(pobj)];%#ok
                            else
                                pobj = reshape(pobj,1,[]);
                                str = [str p.Name ': ' num2str(pobj)];%#ok
                            end
                        elseif isstruct(pobj)
                            if any(size(pobj) > 1)
                                str = [str p.Name ' (struct, [' num2str(size(pobj)) ']), fields: '];%#ok
                            else
                                str = [str p.Name ' (struct), fields: '];%#ok
                            end
                            fn = fieldnames(pobj);
                            for fnidx = 1:length(fn)
                                %str = [str fn{fnidx} this.getObjectConfig(pobj, numtabs+1, depth-1)];%#ok
                                str = [str fn{fnidx} ','];%#ok
                            end
                            str = str(1:end-1);
                        elseif isa(pobj,'function_handle')
                            str = [str p.Name ' (f-handle): ' func2str(pobj)];%#ok
                        end
                    else
                        if isequal(obj,pobj)
                            str = [str p.Name ': [self-reference]'];%#ok
                        else
                            str = [str p.Name ': empty'];%#ok
                        end
                    end
                    str = [str '\n'];%#ok
                    %str = strcat(str,'\n');
                end
            end
            % Take away the last \n
            str = str(1:end-2);
            
            % Format!
            str = sprintf(str);
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
    end
    
end

