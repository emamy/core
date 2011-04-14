classdef Dictionary < handle
    % A basic dictionary of key/value pairs.
    %
    % Usage notes: Despite the standard set/get methods one can also use
    % the matlab-like subscripted notation to set and get values.
    %
    % Example:
    % @code
    % d = general.collections.Dictionary;
    % d('somekey') = 'somevalue';
    % d('nr') = 2346;
    % d('myclass') = KerMorObject;
    % disp(d('nr'));
    % disp(d('myclass'));
    % disp(d('nonexistent'));
    % @endcode
    % 
    % @author Daniel Wirtz @date 2011-04-06
    %
    % @new{0,3,dw,2011-04-06} Added this class. Used e.g. for the
    % customized save/load process to keep integrity of cross-class
    % references (Update 2011-04-07: Customized save/load discontinued, so
    % this class is not used at the moment)
    
    properties(Access=private)
        % The internal List structure
        List;
    end
    
    properties(Dependent)
        %Keys;
        
        %Values;
        
        Count;
    end
    
    methods
        function this = Dictionary
            this.List = struct('Key',{},'Value',{});
        end
        
        function set(this, key, value)
            this.set_(key, value);
        end
        
        function value = get(this, key)
            value = this.get_(key);
        end
        
        function clear(this)
            % Clears the dictionary
            this.List = struct('Key',{},'Value',{});
        end
        
        function value = subsref(this, key)
            % Implements subscripted value retrieval.
            %
            % See also: subsref
            if (isequal(key.type,'()'))
                value = this.get_(key.subs{1});
            else
                builtin('subsref', this, key);
            end
        end
        
        function this = subsasgn(this, key, value)
            % Implements subscripted assignment.
            %
            % See also: subsasgn
            if (isequal(key.type,'()'))
                this.set_(key.subs{1}, value);
            else
                builtin('subsasgn', this, key, value);
            end
        end
        
        function c = get.Count(this)
            c = length(this.List);
        end
    end
        
    methods(Access=private)
        function value = get_(this, key)
            % Default: []
            value = [];
            if isempty(key)
                return;
            elseif ~isa(key, 'char')
                error('The key must be a string');
            end            
            for i=1:length(this.List)
                if isequal(this.List(i).Key,key)
                    % Return found value
                    value = this.List(i).Value;
                    return;
                end
            end
        end
        
        function set_(this, key, value)
            cl = class(value);
            fprintf('Set class %s for key %s\n',cl,key);
            
            for i=1:length(this.List)
                if isequal(this.List(i).Key,key)
                    % Set value if already in list
                    this.List(i).Value = value;
                    return;
                end
            end
            % Otherwise: Extend!
            this.List(end+1).Key = key;
            this.List(end).Value = value;
        end
    end
    
end

