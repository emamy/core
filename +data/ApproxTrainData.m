classdef ApproxTrainData < handle
% ApproxTrainData: Data class for approximation training data, containing
% several useful bounding box properties etc.
%
% @author Daniel Wirtz @date 2011-11-02
%
% @new{0,6,dw,2012-07-18} Moved computation method for approx train data from
% models.BaseFullModel to here.
%
% @new{0,6,dw,2011-11-16} Added a new method data.ApproxTrainData.getCombinedData which returns
% the triples `x_i,t_i,\mu_i` as one column matrix.
%
% @new{0,5,dw,2011-11-02} Added this class.
%
% @todo Make xi,ti,mui private fields and update dependents like Center upon re-setting one of
% those values! (consistency issue)
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % The state space samples `x_i = x(t_i;\mu_i)`, stored row-wise in a data.FileMatrix
        % instance.
        %
        % @type data.FileMatrix
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
        
        % The evaluations of `f(x_i,t_i,\mu_i)`, stored row-wise in a data.FileMatrix
        %
        % @type data.FileMatrix
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
            % Creates a new approx train data instance
            %
            % Parameters:
            % xi: The state space samples `x_i` @type matrix<double>
            % ti: The time instances `t_i` @type rowvec<double>
            % mui: The parameter samples `\mu_i` @type matrix<double>
            %
            % Leave `t_i` and/or `\mu_i` empty if no time or parameter
            % training data is present.
            %
            % Assign local variables
            if ~isa(xi,'data.FileMatrix')
                xi = data.FileMatrix(xi);
            end
            this.xi = xi;
            if isempty(ti)
                ti = double.empty(0,size(xi,2));
            end
            this.ti = ti;
            if isempty(mui)
                mui = double.empty(0,size(xi,2));
            end
            this.mui = mui;
            
            % Build box values
            box = struct;
            [box.xmin, box.xmax] = xi.getColBoundingBox;
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
                [box.mumin, box.mumax] = Utils.getBoundingBox(mui);
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
            % c: A column vector composed as `[x; t; mu]` @type colvec<double>
            % idx: The index in xi,ti,mui samples @type integer
            A = repmat(this.Center, 1, size(this.xi,2));
            B = [this.xi; this.ti; this.mui];
            [~, idx] = min(sum((A-B).^2,1));
            c = B(:,idx);
        end        
       
        function satd = subset(this, sel)
            if isscalar(sel)
                sel = 1:sel:size(this.xi,2);
            end
            xi = this.xi(:,sel);
            ti = [];
            if ~isempty(this.ti)
                ti = this.ti(:,sel);
            end
            mui = [];
            if ~isempty(this.mui)
                mui = this.mui(:,sel);
            end
            satd = data.ApproxTrainData(xi, ti, mui);
            satd.fxi = this.fxi(:,sel);
        end
        
        function [train, validation, randidx] = splitToTrainValidationSets(this, perc, seed)
            % Creates a training and validation data set from this approx
            % train data instance.
            s = RandStream('mt19937ar','Seed',seed);
            randidx = s.randperm(size(this.xi,2));
            valnum = round(size(this.xi,2)*perc);
            train = this.subset(valnum+1:size(this.xi,2));
            validation = this.subset(randidx(1:valnum));
        end
        
        function [e,lbl,pt] = getErrorsFor(this, fun)
            fxi = this.fxi.toMemoryMatrix;
            diff = fxi - fun.evaluate(this.xi.toMemoryMatrix);
            fxino = Norm.L2(this.fxi);
            l2 = Norm.L2(diff);
            lbl{1} = 'Max abs l2';
            e(1) = max(l2);
            lbl{2} = 'Max rel l2';
            e(2) = max(l2./fxino);
            lbl{3} = 'Mean abs l2';
            e(3) = mean(l2);
            lbl{4} = 'Mean rel l2';
            e(4) = mean(l2./fxino);
            
            l1 = Norm.L1(diff);
            fxino = Norm.L1(fxi);
            lbl{5} = 'Max abs l1';
            e(5) = max(l1);
            lbl{6} = 'Max rel l1';
            e(6) = max(l1./fxino);
            lbl{7} = 'Mean abs l1';
            e(7) = mean(l1);
            lbl{8} = 'Mean rel l1';
            e(8) = mean(l1./fxino);
            
            pt = PrintTable;
            pt.Caption = sprintf('Errors on ApproxTrainData (#xi = %d)',size(this.xi,2));
            pt.HasRowHeader = true;
            for k = 1:8
                pt.addRow(lbl{k}, e(k), {'%4.2g'});
            end
            pt = pt';
            
            if nargout < 3
                pt.print;
            end
        end
        
        function relocate(this, new_dir)
            % Convenience method. Relocates the xi and fxi FileMatrix instances if present.
            %
            % Parameters:
            % new_dir: The new directory @type char
            %
            % @new{0,7,dw,2013-05-28} Added this method.
            if ~isempty(this.xi) && isa(this.xi,'data.FileData')
                this.xi.relocate(new_dir);
            end
            if ~isempty(this.fxi) && isa(this.fxi,'data.FileData')
                this.fxi.relocate(new_dir);
            end
        end
        
        function delete(this)
            this.xi = [];
            this.fxi = [];
        end
        
        function makeUniqueXi(this)
            xi = this.xi.toMemoryMatrix;
            [~, selidx] = unique(xi','rows');
            oldnum = size(xi,2);
            if length(selidx) < oldnum
                this.xi = this.xi.spawnWithContent(xi(:,selidx));
                if ~isempty(this.fxi)
                    this.fxi = this.fxi.spawnWithContent(this.fxi(:,selidx));
                end
                if ~isempty(this.ti)
                    this.ti = this.ti(selidx);
                end
                if ~isempty(this.mui)
                    this.mui = this.mui(:,selidx);
                end
                fprintf('makeUniqueXi: Deleted %d training points (old:%d, new:%d)\n',...
                    oldnum-length(selidx),oldnum,length(selidx));
            elseif KerMor.App.Verbose > 0
                fprintf('Nothing to do, all xi unique.\n');
            end
        end
        
        function [minDist, meanDist, maxDist] = getXiDists(this)
            xi = this.xi.toMemoryMatrix;
            xisq = ones(size(xi,2),1)*sum(xi.^2);
            dist = sqrt(xisq' + xisq - 2*(xi'*xi));
            maxDist = max(dist(:));
            meanDist = mean(dist(:));
            dist(logical(eye(size(xi,2)))) = Inf;
            minDist = min(dist(:));
        end
        
        function set.xi(this, value)
            if ~isempty(value)
                if ~isa(value,'data.FileMatrix') && ismatrix(value)
                    value = data.FileMatrix(value);
                elseif ~isa(value,'data.FileMatrix')
                    error('The xi property must be a matrix or a data.FileMatrix');
                end
            end
            this.xi = value;
        end
        
        function set.fxi(this, value)
            if ~isempty(value)
                if ~isa(value,'data.FileMatrix') && ismatrix(value)
                    value = data.FileMatrix(value);
                elseif ~isa(value,'data.FileMatrix')
                    error('The fxi property must be matrix or a data.FileMatrix');
                end
            end
            this.fxi = value;
        end
    end
    
    methods(Static)
        function atd = computeFrom(model, f, selector, parallel)
            % Computes approximation training data for a model, function and selector
            %
            % This method is implemented here as both the original offline
            % phase 4 uses this method but also the DEIM error estimator
            % needs the same routine for the MatrixDEIM approximation.
            %
            % Parameters:
            % model: The full model instance @type models.BaseFullModel
            % f: The model's core function @type dscomponents.ACoreFun
            % selector: The training data selector @type data.selection.ASelector
            % parallel: Flag to set if computation should be done in parallel (MatlabPool)
            % @type logical @default false
            %
            % See also: models.BaseFullModel error.DEIMEstimator
            if nargin < 4
                parallel = false;
            end
            
            % Select subset of projection training data
            atd = selector.selectTrainingData(model);
            
            % If projection is used, train approximating function in
            % centers projected into the subspace.
            projected = false;
            if ~isempty(model.Data.V) && ~isempty(model.Data.W)
                hlp = model.Data.V*(model.Data.W'*atd.xi);
                projected = true;
                % Precautionary case: If the atd.xi data is so small that it fits into one
                % block, the arithmetic operations will return a simple double matrix, loosing
                % and directory information of the original instance. Hence, in this case we
                % create a new instance using the path of the old matrix
                % The goal is to provide filesystem instances inside the model.Data class that
                % are exclusively located inside the model's data folder.
                if ~isa(hlp,'data.FileMatrix')
                    atd.xi = atd.xi.spawnWithContent(hlp);
                else
                    atd.xi = hlp;
                end
            end
            
            % Compute f-Values at training data
            % This only needs to be done if projection is used within the
            % model or the fxi data has not been provided yet by the
            % selector.
            if projected || isempty(atd.fxi)
                if parallel
                    fval = zeros(size(atd.xi));
                    if KerMor.App.Verbose > 0
                        fprintf('Starting parallel f-values computation at %d points on %d workers...\n',size(atd.xi,2),matlabpool('size'));
                    end
                    parfor sidx=1:size(atd.xi,2)
                        fval(:,sidx) = ...
                            f.evaluateMulti(atd.xi(:,sidx),... % x 
                            atd.ti(sidx),... % t
                            atd.mui(:,sidx)); %#ok<PFBNS> % mu
                    end
                    atd.fxi = fval;
                else
                    xi = atd.xi; %#ok<*PROP>
                    if KerMor.App.Verbose > 0
                        fprintf('Serial computation of f-values at %d points (%d xi-blocks) ...\n',size(xi,2),xi.nBlocks);
                    end
                    % Use the same storage location as of the xi FileMatrix
                    fxi = data.FileMatrix(f.fDim,size(xi,2),'Dir',fileparts(xi.DataDirectory));
                    if KerMor.App.Verbose > 1
                        pi = ProcessIndicator('Computing %d %dx%d-block f evaluations on %d %dx%d-blocks of xi snapshots',...
                            fxi.nBlocks,false,fxi.nBlocks,fxi.n,fxi.bCols,xi.nBlocks,xi.n,xi.bCols);
                    end
                    for i=1:fxi.nBlocks
                        pos = fxi.getBlockPos(i);
                        mui = [];
                        if ~isempty(atd.mui)
                            mui = atd.mui(:,pos);
                        end
                        hlp = f.evaluateMulti(xi(:,pos), atd.ti(pos), mui);
                        fxi(:,pos) = hlp;
                        if KerMor.App.Verbose > 1
                            pi.step;
                        end
                    end
                    if KerMor.App.Verbose > 1
                        pi.stop;
                    end
                    atd.fxi = fxi;
                end
            end
        end
        
        function [atd, minScale, maxScale] = scaleXiZeroOne(atd)
            error('not yet implemented with handle usage');
            minScale = min(atd.xi,[],2);
            maxScale = max(atd.xi,[],2);
            n = size(atd.xi,2);
            atd.xi = (atd.xi - repmat(minScale,1,n)) ./ repmat(maxScale - minScale,1,n);
        end
    end
    
end