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
    % @change{0,3,dw,2011-04-20} - Improved the subsref and subsasgn methods to forward eventual
    % further sub-assignments to the respective underlying values. Now i.e. assignments of the type
    % @code d('t').Somefield @endcode natively creates a struct at @code d('t') @endcode.
    % - Added a method general.collections.Dictionary.containsKey
    % - Added two properties general.collections.Dictionary.Keys and
    % general.collections.Dictionary.Values
    %
    % @new{0,3,dw,2011-04-06} Added this class.
    
    properties(Access=private)
        % The internal List structure
        List;
    end
    
    properties(Dependent)
        Keys;
        
        Values;
        
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
        
        function bool = containsKey(this, key)
            if ~isa(key,'char')
                error('Key must be a string.');
            end
            bool = false;
            for i=1:length(this.List)
                if strcmp(this.List(i).Key,key)
                    bool = true;
                    return;
                end
            end
        end
        
        function value = subsref(this, key)
            % Implements subscripted value retrieval.
            %
            % See also: subsref
            if (isequal(key(1).type,'()'))
                value = this.get_(key(1).subs{1});
                if length(key) > 1
                    value = builtin('subsref', value, key(2:end));
                end
            else
                value = builtin('subsref', this, key);
            end
        end
        
        function this = subsasgn(this, key, value)
            % Implements subscripted assignment.
            %
            % See also: subsasgn
            if (isequal(key(1).type,'()'))
                if length(key) > 1
                    value = builtin('subsasgn', this.get_(key(1).subs{1}), key(2:end), value);
                end
                this.set_(key(1).subs{1}, value);
            else
                builtin('subsasgn', this, key, value);
            end
        end
        
        function c = get.Count(this)
            c = length(this.List);
        end
        
        function k = get.Keys(this)
            k = {this.List(:).Key};
        end
        
        function k = get.Values(this)
            k = {this.List(:).Value};
        end
        
        function display(this)
            fprintf('general.collections.Dictionary with %d elements.\n',this.Count);
            arrayfun(@(kv)(disp({'Key:' kv.Key 'Value:' kv.Value})),this.List);
        end
    end
        
    methods(Access=private)
        function value = get_(this, key)
            % Default: []
            value = [];
            if isempty(key)
                return;    
            end
            if isa(key, 'char')
                for i=1:length(this.List)
                    if isequal(this.List(i).Key,key)
                        % Return found value
                        value = this.List(i).Value;
                        break;
                    end
                end
                return;
            elseif isposintscalar(key)
                if key < 1 || key > this.Count
                    error('Key index exceeds dictionary.');
                end
                value = this.List(key).Value;
                return;
            end
            error('The key must be a string or integer');
        end
        
        function set_(this, key, value)
            %cl = class(value);
            %fprintf('Set class %s for key %s\n',cl,key);
            if isa(key, 'char')
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
            elseif isposintscalar(key)
                if key < 1 || key > this.Count
                    error('Setting per numerical index only allowed within existing items; use a char key for automatic insertions.');
                end
                this.List(key).Value = value;
                return;
            end   
        end
    end
    
end

