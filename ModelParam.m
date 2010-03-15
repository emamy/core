classdef ModelParam
    %MODELPARAM Summary of this class goes here
    %   Detailed explanation goes here
    
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
    end
    
    methods
        function this = ModelParam(name, range, desired)
            % Creates a new model parameter.
            %
            % Arguments:
            % range  - can be either a scalar or a 1x2 double vector.
            %
            % @TODO: Validity checks
            
            if nargin > 0
                
                % Double the range in case a scalar is passed
                if isscalar(range)
                    range = [range range];
                end
                
                this.Name = name;
                this.Range = range;
                this.Desired = desired;
            end
        end
        
        function value = get.MinVal(this)
            value = this.Range(1);
        end
        
        function value = get.MaxVal(this)
            value = this.Range(2);
        end
    end
    
end

