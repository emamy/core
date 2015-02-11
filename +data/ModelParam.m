classdef ModelParam < handle
% Stores model parameters.
% 
% This is a simple container for model parameters that stores a
% parameters name and range as well as other parameter-related values.
%
% @author Daniel Wirtz @date 2010-05-01
%
% @new{1,0,dw,2014-09-10} Added spacing to constructor arguments
%
% @change{0,5,dw,2011-09-15} Added some documentation
%
% @change{0,3,sa,2011-05-10} Implemented setters for the properties
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % The Name of the Parameter
        %
        % @default New Parameter @type char
        Name = 'New Parameter';
        
        % The range of the values the parameter may take
        %
        % @default [] @type double
        Range;
        
        % For Sampling: The desired number of samples.
        % This field may be used differently, refer to the sampling module
        % for its usage.
        %
        % @default 1 @type integer
        %
        % See also: sampling
        Desired = 1;
        
        % Sets the desired sampling type for random or grid sampling.
        %
        % Allowed values are 'lin' (linear) and 'log' (logarithmic w.r.t
        % base 10)
        %
        % @type char @default 'lin'
        Spacing = 'lin';
        
        % The default value of this parameter.
        %
        % Must be passed to constructor.
        %
        % @type double
        Default;
    end
    
    properties(Dependent)
        % The maximum value of the parameter's data.ModelParam.Range
        %
        % @type double
        MinVal;
        
        % The minimum value of the parameters Range
        %
        % @type @double
        MaxVal;
        
        % Flag that indicates if this parameter is constant or has a range
        % to vary within.
        %
        % @type boolean
        HasRange;
    end
    
    methods
        function this = ModelParam(name, default, varargin)
            % Creates a new model parameter.
            %
            % Parameters:
            % name: Parameter name. @type char
            % default: The default value @type double
            % varargin: Extra specifications
            % Range: Can be either a scalar or a 1x2 double vector. @type
            % double
            % Desired: The desired number for GridSampling @type integer
            % Spacing: The intended sample spacing over the range @type
            % char @default 'lin'
            %
            % If an argument is specified, all have to be specified. This
            % is only done to enable creation of empty ModelParam-instances
            % for cell arrays, for example.
            
            this.Name = name;
            this.Default = default;
            
            i = inputParser;
            i.addParamValue('Desired',1);
            i.addParamValue('Range',[default * .9 default*1.1]);
            i.addParamValue('Spacing','lin');
            i.parse(varargin{:})
            r = i.Results;            
            this.Desired = r.Desired;
            this.Range = r.Range;
            this.Spacing = r.Spacing;
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
            if value < 0 || ~isscalar(value)
                error('Desired must be a positive integer greater than zero.');
            end
            this.Desired = value;
        end
        
        function set.Default(this, value)
            if isempty(value) || ~isscalar(value)
                error('Default must be a scalar double');
            end
            this.Default = value;
        end
        
        function set.Spacing(this, value)
            if ~ischar(value) || ~any(strcmp(value,{'lin','log'}))
                error('Spacing can either be "lin" or "log".');
            end
            this.Spacing = value;
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

