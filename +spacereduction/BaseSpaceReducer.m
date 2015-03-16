classdef BaseSpaceReducer < KerMorObject & IReductionSummaryPlotProvider
    % Base class for all space reduction algorithms.
    %
    % @author Daniel Wirtz
    % @date 11.03.2010
    
    properties(SetObservable)
        % Flag if to include `f(x_i)` values for each `x_i` value, too
        %
        % @propclass{important} Including this data in the subspace can reduce projection
        % errors of nonlinearities
        %
        % @type logical @default false
        IncludeTrajectoryFxiData = false;
        
        IncludeFiniteDifferences = false;
        
        IncludeBSpan = false;
        
        IncludeAxData = false;
        
        % This property determines which dimensions are to be reduced.
        % It is possible to reduce the entire number of DoFs but also only
        % selected ones.
        %
        % There are three options:
        % - ':' uses all dimensions as target dimensions
        % - A single column index vector determines one subset of DoFs that
        % should be reduced
        % - A cell array of index column vectors gives a list of separate
        % subspaces that should be generated during the state space
        % reduction process. The resulting projection matrices will have a
        % block structure. Overlapping indices are not allowed and will
        % throw an error.
        %
        % @type varargin @default ':'
        TargetDimensions = ':';
    end
    
    properties(SetAccess=protected)
        ProjectionError = [];
    end
    
    methods
        
        function [V, W] = generateReducedSpace(this, model)
            % Generates the reduced linear subspace for the given model.
            %
            % Depending on TargetDimensions, either all DoFs are reduced or
            % one or more subsets, if e.g. physical quantities have to be
            % treated separately.
            %
            % See also: TargetDimensions
            
            % Sanity checks for property are done in setter (below)
            td = this.TargetDimensions;
            % Simplest case: Reduce all DoFs.
            if strcmp(td,':')
                [V,W] = this.generateReducedSpaceImpl(model, ':');
            else
                % For single selections, create a single element cell
                if ~iscell(td)
                    td = {td};
                end
                dim = model.System.StateSpaceDimension;
                V = zeros(dim,0);
                W = zeros(dim,0);
                offset = 0;
                done = false(dim,1);
                e = cell(length(td),1);
                for k=1:length(td)
                    subset = td{k};
                    [Vpart, Wpart] = this.generateReducedSpaceImpl(model, subset);
                    e{k} = this.ProjectionError;
                    extdim = size(Vpart,2);
                    V(subset,offset + (1:extdim)) = Vpart;
                    if ~isempty(W)
                        W(subset,offset + (1:extdim)) = Wpart;
                    end
                    offset = offset + extdim;
                    done(subset) = true;
                end
                this.ProjectionError = e;
                
                % Check if there were dimensions that have not been selected
                % for reduction and augment the projection matrices
                % accordingly (by identity vectors). This is used e.g. if the
                % ODE has algebraic conditions (DAEs)
                V = fill(V);
                W = fill(W);
            end
            
            function Afull = fill(A)
                Afull = A;
                if ~isempty(A) && any(~done)
                    unreduced = 1:dim;
                    unreduced(done) = [];
                    nunred = length(unreduced);
                    
                    basesize = size(A,2);
                    Afull = zeros(dim,basesize+nunred);
                    Afull(done,1:basesize) = A(done,:);
                    % Insert identity for remaining dimensions
                    Afull(~done,basesize+1:basesize+nunred) = eye(nunred);
                end
            end
        end
        
        function plotSummary(this, pm, context)
            if ~isempty(this.ProjectionError)
                e = this.ProjectionError;
                if ~iscell(e)
                    e = {e};
                end
                extra = '';
                for k=1:length(e)
                    if length(e) > 1
                        extra = sprintf(' (Subset %d)',k);
                    end
                    str = sprintf('%s/%s: Projection error over training data%s',...
                    context,class(this),extra);
                    h = pm.nextPlot(sprintf('spacereduction_projerr_%d',k),...
                        str,'subspace size','error');
                    semilogy(h,e{k},'LineWidth',2);
                end
            end
        end
    end
    
    methods(Access=protected)
        
        function V = getInitialSpace(this, blockdata, pod, subset)
            % Computes the initial space, which is the first POD mode of
            % the initial values!
            
            n = blockdata.getNumBlocks;
            x = blockdata.getBlock(1);
            x0 = x(subset,1);
            if this.ComputeParallel
                parfor idx=2:n
                    x = blockdata.getBlock(idx);%#ok
                    x0 = [x0, x(subset,1)];
                end
                x0 = unique(x0','rows')';
            else
                for idx=2:n
                    x = blockdata.getBlock(idx);
                    x = x(subset,1);
                    % Only add nonexisting vectors
                    if isempty(Utils.findVecInMatrix(x0,x))
                        x0 = [x0 x];%#ok
                    end
                end
            end
            if all(x0 == 0)
                if KerMor.App.Verbose > 1
                    fprintf('Initial values are all zero vectors. Using main POD mode of first block data as initial space.\n');
                end
                xb1 = blockdata.getBlock(1);
                V = pod.computePOD(xb1(subset,:));
            elseif size(x0,2) > 1
                V = pod.computePOD(x0);
            else
                V = x0;
            end
            V = V / norm(V);
        end
    end
        
    methods(Access=protected, Abstract)
        [V,W] = generateReducedSpaceImpl(this, model, subset);
    end
    
    %% Getter & setter
    methods 
        function set.IncludeTrajectoryFxiData(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeTrajectoryFxiData = value;
        end
        
        function set.IncludeFiniteDifferences(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeFiniteDifferences = value;
        end
        
        function set.IncludeBSpan(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeBSpan = value;
        end
        
        function set.IncludeAxData(this, value)
            if isempty(value) || ~islogical(value) || ~isscalar(value)
                error('Value needs to be a logical scalar.');
            end
            this.IncludeAxData = value;
        end
        
        function set.TargetDimensions(this, value)
            err = false;
            if ~strcmp(value ,':')
                if iscell(value)
                    for k=1:length(value)
                        if ~isvector(value{k})
                            err = true;
                            break;
                        end
                    end
                elseif ~isvector(value)
                    err = true;
                end
            end
            if err
                error('TargetDimensions must either be ":", a column vector or a cell of column vectors.');
            end
            this.TargetDimensions = value;
        end
    end
    
end

