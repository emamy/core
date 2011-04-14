classdef NNTracker < handle
    % Nearest Neighbor tracking class
    %
    % This class allows to successively build a list of nearest neighbor
    % distances for a set of points. Call the general.NNTracker.addPoint
    % method to insert a new point. This automatically updates all existing
    % points' distances.
    % 
    % @author Daniel Wirtz @date 2011-03-30
    % 
    % @new{0,3,dw,2011-03-30} Added this class to support nearest neighbor
    % tracking (in `L^2`-norm) for a successively built set of points.
    %
    % @new{0,3,dw,2011-04-01} Added a UniqueValues property that determines
    % if repeatedly added values are ignored or not.
    
    properties(SetAccess=private)
        % The points matrix containing each point as a column vector
        % For each point `x_i` with column index `i` the corresponding
        % nearest neighbor distance `NN(x_i)` can be found in the property
        % NNDists(i)
        %
        % See also: NNDists
        Points;
        
        % The nearest neighbor distances `NN(x_i)` for each point `x_i`
        % from the property Points.
        %
        % See also: Points
        NNDists;
        
        % Flag that determines if values that are added more than once are
        % ignored or not.
        %
        % @default true
        UniqueValues = true;
    end
    
    methods
        
        function this = NNTracker
            % Creates a new NNTracker instance
            this.Points = [];
            this.NNDists = [];
        end
        
        function addPoint(this, x)
            % Adds a point to the list of points the nearest neighbors are
            % kept track of.
            %
            % Parameters:
            % x: The point to add as a column vector.
            
            % Checks
            if size(x,2) ~= 1
                error('The point to add must be a column vector.');
            end
            
            np = size(this.Points,2);
            % Previous points exist
            if np > 0
                
                % Prepare comparison matrix
                X = repmat(x,1,np);
                d = sqrt(sum((this.Points-X).^2,1));
                
                % Check if unique values should be taken
                if this.UniqueValues && any(d == 0)
                    if KerMor.App.Verbose > 3
                        fprintf('NNTracker: Double value %s detected\n',num2str(x'));
                    end
                    return;
                end
                
                % Update any closer distances
                sm = d < this.NNDists;
                this.NNDists(sm) = d(sm);
                
                % Add point itself
                this.NNDists(end+1) = min(d);
                
            % First point added
            else
                this.NNDists = Inf;
            end
            this.Points(:,end+1) = x;
        end
        
        function d = getMaxNN(this)
            % Returns the maximum nearest neighbor distance for all
            % currently registered points.
            %
            % Return values:
            % d: The maximum nearest neighbor distance
            d = max(this.NNDists);
        end
        
        function d = getMinNN(this)
            % Returns the maximum nearest neighbor distance for all
            % currently registered points.
            %
            % Return values:
            % d: The maximum nearest neighbor distance
            d = min(this.NNDists);
        end
    end
    
end

