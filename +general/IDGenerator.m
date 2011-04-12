classdef IDGenerator < handle
    % Generates unique IDs
    %
    % Used to assign unique IDs to any KerMorObject instance
    %
    % See also: KerMorObject.ID
    %
    % @author Daniel Wirtz @date 2011-04-06
    %
    % @new{0,3,dw,2011-04-6} Added this class. Used to assign unique IDs to
    % any KerMorObjects.

    properties(SetAccess=private)
        % The seed used as base for new IDs
        Seed;
        
        % Counts how many IDs have been generated since singleton
        % instantiation
        Counter;
    end
    
    methods(Access=private)
        function this = IDGenerator
            % Creates a new IDGenerator.
            %
            % Private constructor as this class is a singleton during its
            % lifetime.
            this.Seed = sum(clock)+cputime;
            this.Counter = 1;
        end
        
        function id = genID(this)
            % @todo Use CalcMD5 here once it compiles
            id = num2str(this.Seed+this.Counter);
            this.Counter = this.Counter+1;
        end
    end
    
    methods(Static)
        function gen = generateID
            % Generates a new unique ID
            %
            % Return values:
            % gen: A string ID
            persistent g;
            if isempty(g)
                g = general.IDGenerator;
            end
            gen = g.genID;
        end
    end
    
end

