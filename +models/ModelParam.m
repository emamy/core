classdef ModelParam < handle
    %MODELPARAM Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @change{0,3,sa,2011-05-10} Implemented setters for the properties
    %
    
    properties
        % The Name of the Parameter
        Name = 'New Parameter';
        
        % The range of the values the parameter may take
        Range;
        
        % For Sampling: The desired number of samples.
        % This field may be used differently, refer to the sampling module
        % for its usage.
        %
        % See also: sampling
        Desired = 1;
    end
    
    properties(Dependent)
        MinVal;
        MaxVal;
        HasRange;
    end
    
    methods
        function this = ModelParam(name, range, desired)
            % Creates a new model parameter.
            %
            % Paramters:
            % name: Parameter name
            % range: Can be either a scalar or a 1x2 double vector.
            % desired: The desired number for GridSampling
            %
            % If an argument is specified, all have to be specified. This
            % is only done to enable creation of empty ModelParam-instances
            % for cell arrays, for example.
            %
            % @todo: Validity checks
            
            if nargin > 0
                this.Name = name;
                this.Range = range;
                this.Desired = desired;
            end
        end
        
        function set.Name(this, value)
            if ~ischar(value)
                error('name should be a character field');
            end
            this.Name = value;
        end
        
        function set.Range(this, range)
            % Double the range in case a scalar is passed
            if isscalar(range)
                range = [range range];
            end
            if range(2) < range(1)
                error('Invalid range: MinVal must be greater or equal to MaxVal.');
            end
            this.Range = range;
        end
        
        function set.Desired(this, value)
            if value < 1 || ~isscalar(value)
                error('Desired must be a positive integer greater than zero.');
            end
            this.Desired = value;
        end     
        
        function value = get.MinVal(this)
            value = this.Range(1);
        end
        
        function value = get.MaxVal(this)
            value = this.Range(2);
        end
        
        function value = get.HasRange(this)
            value = abs(this.MinVal - this.MaxVal) > 10*eps;
        end
    end
    
end

