classdef DefaultSelector < approx.selection.ASelector
% Selects all training data.
% 
% @author Daniel Wirtz @date 2011-04-12
%
% @new{0,4,dw,2011-05-06}
% - Implemented ICloneable interface.
%
% @new{0,3,dw,2011-04-12} Added this default selection algorithm.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    methods
        function copy = clone(this)
            copy = approx.selection.DefaultSelector;
            copy = clone@approx.selection.ASelector(this, copy);
        end
    end

    methods(Access=protected,Sealed)
        function [xi, ti, mui] = select(this, model)%#ok
            xi = [];
            ti = [];
            mui = [];
            for k=1:model.Data.getNumTrajectories
                [x, mu] = model.Data.getTrajectoryNr(k);
                xi = [xi x]; %#ok
                ti = [ti model.Times]; %#ok
                mui = [mui repmat(mu,1,size(x,2))]; %#ok
            end
        end
    end
    
end