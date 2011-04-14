classdef LinspaceSelector < approx.selection.ASelector
% Selects Size equally spaced samples from the training data.
% 
% Equally spaced in this context means with respect to the indices, NOT to the spatial distances
% between the training samples.
%
% @author Daniel Wirtz @date 2011-04-12
%
% @new{0,3,dw,2011-04-12} Added this class to allow for homogeneous training data selection.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The (maximum) number of elements to take
        %
        % @default `\infty`
        Size = Inf;
    end
    
    methods(Access=protected,Sealed)
        function atd = select(this, model)
            % Selects Size equally spaced samples from the training data.
            % 
            % Equally spaced in this context means with respect to the
            % indices, NOT to the spatial distances between the training
            % samples.
            sn = model.Data.TrainingData;
            s = min(this.Size,size(sn,2));
            selection = round(linspace(1,size(sn,2),s));
            atd = sn(:,selection);
            this.LastUsed = selection;
        end
    end
    
    methods
        function set.Size(this, value)
            if ~isposintscalar(value)
                error('The value must be a positive integer.');
            end
            this.Size = value;
        end
    end
    
end

