classdef ApproxTrainData
% ApproxTrainData: Data class for approximation training data, containing
% several useful bounding box properties etc.
%
% @author Daniel Wirtz @date 2011-11-02
%
% @new{0,6,dw,2011-11-16} Added a new method data.ApproxTrainData.getCombinedData which returns
% the triples `x_i,t_i,\mu_i` as one column matrix.
%
% @new{0,5,dw,2011-11-02} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The state space samples `x_i = x(t_i;\mu_i)`
        %
        % @type matrix
        xi = [];
        
        % The time samples `t_i`
        %
        % @type rowvec
        ti = [];
        
        % The parameter samples `\mu_i` used computing the parent
        % trajectories of `x_i`
        %
        % @default [] @type matrix
        mui = [];
        
        % The evaluations of `f(x_i,t_i,\mu_i)`
        %
        % @type matrix
        fxi;
    end
    
    properties(SetAccess=private)
        % The geometrical center of the datas bounding box
        %
        % A column vector composed as `[x; t; mu]`
        % @type colvec
        Center;
        
        % The diameter of the state space samples
        % @type double
        xiDia;
        
        % The time span of the time samples
        % @type double
        tiDia;
        
        % The diameter of the parameter space samples
        % @type double
        muiDia;
        
        % Flag that indicates if time samples are present
        % @type logical @default false
        hasTime = false;
        
        % Flag that indicates if param samples are present
        % @type logical @default false
        hasParams = false;
        
        % The index of the time entry in the combined `[x;t;\mu]` vector.
        %
        % Zero means no time samples are present.
        %
        % @type integer @default 0
        tOff = [];
        
        % The index of the first parameter entry in the combined
        % `[x;t;\mu]` vector.
        %
        % Zero means no parameter samples are present.
        %
        % @type integer @default 0
        muOff = [];
        
        % A box struct allowing access to the specific min and max values
        % of the samples. Recommended for Debug use only.
        %
        % @type struct
        Box = [];
    end
    
    methods
        function this = ApproxTrainData(xi, ti, mui)
            % Assign local variables
            this.xi = xi;
            this.ti = ti;
            this.mui = mui;
            
            % Build box values
            box = struct;
            [box.xmin, box.xmax] = general.Utils.getBoundingBox(xi);
            this.Center = (box.xmin+box.xmax)/2;
            % Get bounding box diameters
            this.xiDia = norm(box.xmax - box.xmin);
            if ~isempty(ti)
                this.hasTime = true;
                box.tmin = min(ti); box.tmax = max(ti);
                this.tiDia = box.tmax-box.tmin;
                this.Center = [this.Center; (box.tmin+box.tmax)/2];
                this.tOff = size(xi,1)+1;
            end   
            if ~isempty(mui)
                this.hasParams = true;
                [box.mumin, box.mumax] = general.Utils.getBoundingBox(mui);
                this.muiDia = norm(box.mumax - box.mumin);

                this.Center = [this.Center; (box.mumin + box.mumax)/2];
                this.muOff = size(xi,1)+1;
                if ~isempty(this.tOff)
                    this.muOff = this.muOff + 1;
                end
            end
            this.Box = box;
        end
        
        function [c, idx] = getClosestToCenter(this)
            % The sample triple closest to the Center
            %
            % Return values:
            % c: A column vector composed as `[x; t; mu]` @type colvec
            % idx: The index in xi,ti,mui samples @type integer
            A = repmat(this.Center, 1, size(this.xi,2));
            B = [this.xi; this.ti; this.mui];
            [~, idx] = min(sum((A-B).^2,1));
            c = B(:,idx);
        end
        
        function comb = getCombinedData(this)
            % Returns the training data as matrix, composed as '[xi; ti; mui]'.
            %
            % Return values:
            % comb: The training data matrix. @type matrix
            comb = this.xi;
            if this.hasTime
                comb = [comb; this.ti];
            end
            if this.hasParams
                comb = [comb; this.mui];
            end
        end
    end
    
end