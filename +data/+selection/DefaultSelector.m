classdef DefaultSelector < data.selection.ASelector
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    methods
        function copy = clone(this)%#ok
            copy = data.selection.DefaultSelector;
            %copy = clone@data.selection.ASelector(this, copy);
        end
    end

    methods(Access=protected,Sealed)
        function [xi, ti, mui] = select(this, model)%#ok
            % Selects ALL the trajectory data as training.
            %
            % Parameters:
            % model: The full model with the training data @type models.BaseFullModel
            %
            % Return values:
            % xi: The selected `x_i = x(t_i)` training data @type data.FileMatrix
            % ti: The selected training times `t_i` @type rowvec
            % mui: The selected parameter samples `\mu_i` with which the states
            % `x_i` have been reached @type matrix
            
            md = model.Data;
            td = md.TrajectoryData;
            nt = td.getNumTrajectories;
            xi = [];
            ti = [];
            mui = [];
            len = length(model.Times);
            if nt > 0
                [xdim, mudim] = td.getTrajectoryDoFs;
                % Use 512 MB chunks for approx train data
                xi = data.FileMatrix(xdim,nt*len,'Dir',md.DataDirectory,'BlockSize',512);
                ti = zeros(1,nt*len);
                if mudim > 0
                    mui = zeros(mudim,nt*len);
                end
                for k=1:nt
                    [x, mu] = td.getTrajectoryNr(k);
                    pos = (k-1)*len+1:k*len;
                    xi(:,pos) = x;
                    ti(pos) = model.Times;
                    if mudim > 0
                        mui(:,pos) = repmat(mu,1,size(x,2));%#ok
                    end
                end
            end
        end
    end
    
end