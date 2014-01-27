classdef Domain < handle
% Domain: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-01-15
%
% @new{0,7,dw,2014-01-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % Samples from the full parameter domain box specified by the
        % parameter's MinVal and MaxVal.
        % Points(:,Locations) specifies the domain in which the parameter
        % samples must lie.
        %
        % Dimension: `d_p \times n_s`, where `n_s` is the number of samples
        % and `d_p` is the dimension of the parameter space, i.e. the
        % system's ParamCount.
        %
        % @type matrix<double> @default []
        Points;
        
        % Logical vector specifying which of the points in Points
        % are valid.
        %
        % Must have the same number of columns as Points.
        %
        % @type vector<logical> @default []
        Locations;
    end
    
    methods
        function this = Domain(points, locations)
            if size(points,2) ~= size(locations,2)
                error('Column count of point matrix and number of location flags must be the same.');
            end
            this.Points = points;
            this.Locations = locations;
        end
        
        function [params, idx] = filter(this, params)
            % Filters from the given parameters only those that belong to
            % this domain
            %
            % Parameters:
            % params: The parameters to filter @type matrix<double>
            %
            % Return values:
            % params: The subset of params which lie inside the Domain
            % @type matrix<double>
            % idx: The indices of the chosen parameters within params @type
            % rowvec<integer>
            nn = dsearchn(this.Points', params');
            idx = find(this.Locations(nn));
            params = params(:, idx);
        end
        
        function set.Locations(this, value)
            if ~islogical(value)
                error('Locations must be a logical row vector.');
            end
            this.Locations = value;
        end
    end
    
end