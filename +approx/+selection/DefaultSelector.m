classdef DefaultSelector < approx.selection.ASelector
% Selects all training data.
% 
% @author Daniel Wirtz @date 2011-04-12
%
% @new{0,3,dw,2011-04-12} Added this default selection algorithm.
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    methods(Access=protected,Sealed)
        function atd = select(this, model)
            atd = model.Data.TrainingData;
            this.LastUsed = 1:size(atd,2);
        end
    end
    
end