classdef ManualSampler < sampling.BaseSampler
% ManualSampler: Allows to set parameter samples for the reduction process
%
% @author Daniel Wirtz @date 2012-04-26
%
% @new{0,6,dw,2012-04-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The samples to use for offline generations.
        %
        % @propclass{critical} If manual sampling is used, this property
        % builds the basic set of parameter samples to be used for offline
        % computations and hence must be chosen with maximum care.
        %
        % @type matrix<double> @default []
        Samples = [];
    end
    
    methods
        function this = ManualSampler(samples)
            this = this@sampling.BaseSampler;
            this.registerProps('Samples');
            if nargin == 1
                this.Samples = samples;
            end
        end
        
        function samples = performSampling(this, ~)
            samples = this.Samples;
        end
    end
    
end