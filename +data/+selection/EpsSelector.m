classdef EpsSelector < data.selection.ASelector
% EpsSelector: Selects as many points from the data that any trajectory point lies within an epsilon
% radius of a training point.
%
% This selection strategy is applied in the TPWLApprox-Setting, but may also be used elsewhere.
%
% @author Daniel Wirtz @date 2011-05-05
%
% @change{0,5,dw,2011-08-04} Adopted this selector class to the new data.ATrajectoryData structure. For
% now it only works on models that sample a single trajectory, as originally proposed by the TPWL
% guys.
%
% @new{0,4,dw,2011-05-06} Implemented the ICloneable interface.
%
% @new{0,4,dw,2011-05-05} Added this class and included it into the DPCS.
%
% @new{0,4,dw,2011-05-04} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The distance within there has to be an expansion point for each
        % x. Gets multiplied by `\sqrt{d}`, where `d` denotes the spatial
        % dimension of the snapshots (projection training data)
        %
        % @propclass{critical} Determines how many training samples are taken. Value MUST be set
        % taking into account the full system's dimension or at least the bounding box of the samples!
        EpsRad = 5;
        
        % @propclass{experimental}
        SubspaceProject = true;
    end
    
    methods
        function this = EpsSelector
            this = this@data.selection.ASelector;
            this.registerProps('EpsRad','SubspaceProject');
            this.EpsRad = 1;
        end
        
        function copy = clone(this)
            copy = data.selection.EpsSelector;
            %copy = clone@data.selection.ASelector(this, copy);
            copy.EpsRad = this.EpsRad;
            copy.SubspaceProject = this.SubspaceProject;
        end
        
        function set.EpsRad(this, value)
            if ~isposrealscalar(value)
                error('EpsRad must be a positive real scalar.');
            end
            this.EpsRad = value;
        end
    end
    
    methods(Access=protected)
        function [xi, ti, mui] = select(this, model)
            % Selects training points with distance EpsSelector.EpsRad,
            % starting from the initial values for each trajectory.
            %
            % Parameters:
            % model: The full model with the training data @type models.BaseFullModel
            %
            % Return values:
            % xi: The selected `x_i = x(t_i)` training data @type matrix
            % ti: The selected training times `t_i` @type rowvec
            % mui: The selected parameter samples `\mu_i` with which the states
            % `x_i` have been reached @type matrix
            td = model.Data.TrajectoryData;
            if td.getNumTrajectories > 1
                warning('KerMor:approx:selection:EpsSelector','The epsilon-selector can only be used with one trajectory at the moment. Ignoring all but the first.');
            end
            [x, mu] = td.getTrajectoryNr(1); 
            if this.SubspaceProject && ~isempty(model.SpaceReducer)
                W = model.Data.W;
                if isempty(W)
                    W = model.Data.V;
                end
                x = model.Data.V*(W'*x);
            end
            %[m,M] = general.Utils.getBoundingBox(x);
            %d = norm(M-m) / 10;
            %this.EpsRad = d;
            selidx = 1;
            idx = 1; 
            cur = x(:,idx);
            while idx < size(x,2)
                while(idx < size(x,2) && norm(x(:,idx)-cur) < this.EpsRad)
                    idx = idx+1;
                end
                selidx(end+1) = idx;%#ok
                cur = x(:,idx);
                idx = idx+1;
            end
            xi = x(:,selidx);
            ti = model.Times(selidx);
            mui = repmat(mu,1,length(selidx));
        end
    end
    
end