function [t, matches, parents] = FindInstance(obj, type, varargin)
% FindInstance: 
%
%
%
% @author Daniel Wirtz @date 2012-06-11
%
% @new{0,6,dw,2012-06-11} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    if ~isa(obj, 'handle')
        error('boo');
    end
    
    matches = {};
    visited = {};
    parents = {};
    t = PrintTable('Instances of "%s" in %s (%s)',...
        type,inputname(1),class(obj));
    t.HasHeader = true;
    t.addRow('Location','Class','ID','Parent',varargin{:});
    
    warning('off','MATLAB:structOnObject');
    findinst_recur(obj, inputname(1), true);
    warning('on','MATLAB:structOnObject');
    
    function findinst_recur(obj, lvl, checkloop)
        if checkloop && any(cellfun(@(o)isequal(o,obj),visited))
            return;
        end
        
        if numel(obj) > 1
            if iscell(obj)
                for l = 1:length(obj)
                    findinst_recur(obj{l}, sprintf('%s{%d}',lvl,l), false);
                end
            else
                for l = 1:length(obj)
                    findinst_recur(obj(l), sprintf('%s(%d)',lvl,l), false);
                end
            end
        else
            if isa(obj,'handle')
                os = struct(obj);
            elseif isstruct(obj)
                os = obj;
            end
            fn = fieldnames(os);
            for idx = 1:length(fn)
                if ~isempty(os.(fn{idx}))
                    prop = os.(fn{idx});
                    thislvl = [lvl '.' fn{idx}];
                    if isa(prop,type)
                        matches{end+1} = prop;%#ok
                        parents{end+1} = obj;%#ok

                        st = 'none';
                        if any(cellfun(@(o)isequal(o,'ID'),fieldnames(prop)))
                            st = prop.ID;
                        end
                        props = cell(1,length(varargin));
                        for k = 1:length(varargin)
                            props{k} = prop.(varargin{k});
                        end
                        t.addRow(thislvl,class(prop),st,class(obj),props{:});
                    end    
                    if isa(prop, 'handle') || isstruct(prop)
                        visited{end+1} = obj;%#ok
                        findinst_recur(prop, thislvl, true);
                    end
                end
            end
        end
    end
    
    
end