classdef EpsSelector < approx.selection.ASelector
% EpsSelector: Selects as many points from the data that any trajectory point lies within an epsilon
% radius of a training point.
%
% This selection strategy is applied in the TPWLApprox-Setting, but may also be used elsewhere.
%
% @author Daniel Wirtz @date 2011-05-05
%
% @new{0,4,dw,2011-05-06} 
% - Implemented the ICloneable interface.
% @new{0,4,dw,2011-05-05} 
% - Added this class and included it into the DPCS.
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
            this = this@approx.selection.ASelector;
            this.registerProps('EpsRad','SubspaceProject');
            this.EpsRad = 1;
        end
        
        function copy = clone(this)
            copy = approx.selection.EpsSelector;
            copy = clone@approx.selection.ASelector(this, copy);
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
        function atd = select(this, model)
            sn = model.Data.TrainingData;
            x = sn(4:end,:);
            if this.SubspaceProject && ~isempty(model.SpaceReducer)
                x = model.Data.V*(model.Data.W'*x);
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
            atd = sn(:,selidx);
            this.LastUsed = selidx;
        end
    end
    
end