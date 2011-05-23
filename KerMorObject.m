classdef KerMorObject < handle
    % Base class for any KerMor class
    %
    % Every KerMor object is identifiable by a unique ID.
    %
    % Moreover, a messaging system is about to be implemented in KerMor
    % to obtain different modes of output and logging. Every class should
    % have access to that system, therefore a protected method will be
    % available within any subclass. (Current skeleton: logMessage)
    % 
    % @author Daniel Wirtz @date 2011-04-05
    %
    % @new{0,3,dw,2011-04-21} Implemented the property default changed supervision system as
    % described in @ref propclasses. Those changes affect most KerMor classes.
    %
    % @change{0,3,dw,2011-04-07} 
    % - Moved the checkType method from models.BaseModel to this class.
    % - Removed the ALoadable class and moved the ID property here.
    %
    % @new{0,3,dw,2011-04-05} Added this class to give any object a common
    % base class and some shared functionality. For this version it will be
    % loading and message logging, possibly cloning as well.
    
    properties(SetAccess=private)
        % An ID that allows to uniquely identify this KerMorObject
        ID = [];
        
        PropertiesChanged = [];
    end
    
    methods
        function this = KerMorObject
            % Constructs a new KerMor object.
            %
            % Important notice at this stage: Due to possible multiple inheritance any object's
            % constructor should check if any custom properties that are assigned during
            % construction are already present / different from the default value given at property
            % declaration. If so, chances are that the constructor is called again for the same
            % object; in this case, assigning a new ID and property changed dictionary caused any
            % old registered properties to be overwritten!
            
            % Check if a constructor for this object has already been called!
            if isempty(this.ID)
                this.ID = general.IDGenerator.generateID;
                this.PropertiesChanged = general.collections.Dictionary;
            end
        end
    end
    
    methods(Access=protected)
        function registerProps(this, varargin)
            % @docupdate Here!
            %
            % @todo - Check if it makes sense to "link" properties with other known properties that
            % helps to emphasize dependencies. maybe even a callback that validates a connected pair?
            % - write sections for different property classes in documentation
            % - include a disable propchange listening flag for use during simulations.
            % - maybe move the Text and Short property extractions to the printPropertyChangedReport
            % method? -> speedup
            % - include DefaultConfirmed flag into output!
            
            %% Validity checks
            % Find property in class
            if isempty(varargin)
                error('Minimum one argument is required.');
            end
            mc = metaclass(this);
            % Iterate over all properties that are to be registered
            for pidx = 1:length(varargin)
                if ~ischar(varargin{pidx})
                    error('All input arguments must be strings.');
                end
                found = false;
                for idx = 1:length(mc.Properties)
                    if strcmp(mc.Properties{idx}.Name,varargin{pidx})
                        p = mc.Properties{idx};
                        found = true;
                        break;
                    end
                end
                % Validity checks
                if ~found
                    error('Property %s does not exists in class %s.',varargin{pidx},mc.Name);
                elseif ~strcmp(p.GetAccess,'public')
                    error('Property %s has GetAccess=%s, ''public'' is required for registering.',p.Name,p.GetAccess);
                elseif p.Constant
                    error('Property %s is marked Constant.',p.Name,p.GetAccess);
                elseif p.Transient
                    error('Transient property %s is not admissible.',p.Name,p.GetAccess);
                elseif strcmp(p.SetAccess,'private')
                    error('Property %s has SetAccess=private, at least ''protected'' is required for registering.',p.Name,p.GetAccess);
                end
                
                key = [p.DefiningClass.Name '.' p.Name];
                hlp = strtrim(help(key)); %key matches the help search path! see doc help
                n = regexp(hlp, '@propclass{(?<level>\w*)}\s*(?<text>[^@]*)','names');
                % Check if property has @propclass tag set
                if ~isempty(n)
                    % Check validity
                    if ~any(strcmp(n.level,KerMorObject.getPropClasses))
                        error('Invalid property class: %s',n.level);
                    end
                    % Add the listener
                    if ~any(strcmp(n.level,'data'))
                        if ~p.SetObservable
                            error('Non-passive registered property %s must have the SetObservable flag.',p.Name);
                        end
                        addlistener(this,p.Name,'PostSet',@(src,evd)this.PropPostSetCallback(src,evd));
                    end
                    % Only overwrite if not present -> save/load process!
                    if ~this.PropertiesChanged.containsKey(key)
                        % Initialize Changed flag to false
                        ps.Changed = false;
                        
                        % Get current value and save as default
                        ps.Default = this.(p.Name);
                        ps.Name = p.Name;
                        
                        % Add misc help texts and classification
                        ps.Short = general.Utils.getHelpShort(hlp);
                        ps.Level = n.level;
                        ps.Text = strtrim(n.text);
                        
                        this.PropertiesChanged(key) = ps;
                    end
                else
                    error('When registering a property you must define the @propclass{<level>} tag in the properties help text.');
                end
            end
        end
    end
        
    methods(Sealed, Access=protected)
        function checkType(this, obj, type)%#ok
            % Object typechecker.
            %
            % Checks if a given object is of the specified type and throws
            % an error if not. Convenience method that can bu used in subclasses at any setter
            % method that ensures the correct value type.
            if ~isempty(obj) && ~isa(obj, type)
                error(['Wrong type ''' class(obj) ''' for this property. Has to be a ' type]);
            end
        end
    end
    
    methods(Access=private)
        function PropPostSetCallback(this, src, evd)
            p = evd.Source;
            key = [p.DefiningClass.Name '.' p.Name];
            ps = this.PropertiesChanged(key);
            if isempty(ps)
                warning('KerMor:devel','PostSet called on property %s but dict does not contain key',key);
            else
                if ~ps.Changed
                    this.PropertiesChanged(key).Changed = true;
                    this.PropertiesChanged(key).DefaultConfirmed = isequal(ps.Default,evd.AffectedObject.(p.Name));
                    % Save some space!
                    this.PropertiesChanged(key).Default = [];
                    this.PropertiesChanged(key).Text = [];
                end
            end
        end
    end
    
    methods(Static,Access=protected)
        function cl = getPropClasses
            cl = {'experimental','critical','important','alglimit','scaling','optional','data'};
        end
        
        function obj = loadobj(obj)
            % Re-register any registered change listeners!
            keys = obj.PropertiesChanged.Keys;
            for idx = 1:obj.PropertiesChanged.Count
                ps = obj.PropertiesChanged(keys{idx});
                if ~ps.Changed && ~any(strcmp(ps.Level,'data'))
                    addlistener(obj,ps.Name,'PostSet',@(src,evd)obj.PropPostSetCallback(src,evd));
                end
            end
        end
    end

end

