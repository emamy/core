classdef AAdaptiveBase < approx.algorithms.ABase
% Base class for adaptive component-wise kernel approximation algorithms
%
% @author Daniel Wirtz @date 2011-11-02
%
% @new{0,6,dw,2012-01-26} Added a new approximation stallment detection. The properties
% MinImprovePerc and ImproveRange control at which stage the approximation is to be stopped if
% no sufficient progress in approximation error is made.
%
% @new{0,5,dw,2011-11-02} 
% - Created this class. Collects common properties of the adaptive approx algorithms and
% provides convenience methods for subclasses.
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
%
% See also: BaseApprox KernelApprox approx.algorithms.ABase
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    properties(SetObservable)
        % The maximum size of the expansion to produce.
        %
        % Equals the maximum number of iterations to perform during
        % adaptive approximation computation as each iteration yields a new
        % center.
        %
        % @propclass{alglimit} 
        % Some text describing the importance of this property.
        %
        % @type integer @default 300
        MaxExpansionSize = 300;
        
        % Stopping condition property. Maximum relative error that may occur
        %
        % @propclass{critical}
        %
        % @default 1e-5 @type double
        MaxRelErr = 1e-5;
                       
        % Determines how the initial center(s) are chosen.
        %
        % Possible values:
        % - 'maxfx' The center with the largest (w.r.t. the ABase.ErrorFun
        % norm) associated `f(x_i)` value is chosen.
        % - 'center' The training point closest to the geometrical center
        % if the training data's bounding box is used.
        % - 't0' Tries to find a training sample for time zero and uses
        % that. If None are found, the strategy falls back to 'center'. If
        % multiple are found, the one closest to the geometrical center is
        % chosen.
        %
        % @propclass{optional} The first center to choose should not greatly influence the
        % algorithm outcome. The default 'maxfx' is related the closest to the greedy-type
        % approach.
        %
        % @type char @default 'maxfx'
        InitialCenter = 'maxfx';
    end
    
    properties(SetAccess=protected)        
        % @default [] @type colvec<double>
        ExpansionSizes;
        
        % The indices of the effectively used centers during the
        Used = [];
    end
    
    properties(Access=protected)
        initialidx;
    end
    
    methods    
        function this = AAdaptiveBase
            this = this@approx.algorithms.ABase;
            
            % Register default property changed listeners
            this.registerProps('MaxExpansionSize','MaxRelErr','InitialCenter');
        end
                        
        function copy = clone(this, copy)
            % Clones the instance.
            copy = clone@approx.algorithms.ABase(this, copy);
            
            % copy local props
            copy.MaxExpansionSize = this.MaxExpansionSize;
            copy.MaxRelErr = this.MaxRelErr;
            copy.initialidx = this.initialidx;
            copy.InitialCenter = this.InitialCenter;
            copy.Used = this.Used;
        end
    end
    
    methods(Access=protected, Sealed)
        function kexp = templateComputeApproximation(this, atd, atd_val)
            % Performs adaptive approximation generation.
            %
            % Parameters:
            % atd: The approximation training data instance @type
            % data.ApproxTrainData
            % atd_val: A validation data set. Not used with VKOGA yet.
            %
            % Return values:
            % kexp: The kernel expansion approximation. @type
            % kernels.KernelExpansion
            
            %% Checks
            if size(atd.xi,2) < this.MaxExpansionSize
                warning('AAdaptiveBase:expansionsize',...
                    'Only %d training samples but having MaxExpansionSize=%d, changing to %d.',...
                    size(atd.xi,2),this.MaxExpansionSize,size(atd.xi,2));
                this.MaxExpansionSize = size(atd.xi,2);
            end
            
            % Init debug fields
            nc = this.ExpConfig.getNumConfigurations;
            this.MaxErrors = zeros(nc,this.MaxExpansionSize);
            this.MaxRelErrors = this.MaxErrors;
            this.ExpansionSizes = zeros(nc,1);
            
            % Compute initial center
            [~, this.initialidx] = this.getInitialCenter(atd);
            
            % Start adaptive extension part of subclass
            kexp = this.startAdaptiveExtension(atd);
        end
    end
    
    %% Convenience methods for use/override in concrete algorithms
    methods(Access=protected)
        
        function initExpansion(this, kexp, atd)
            % Initializes the Find and set first expansion center
            kexp.clear;
            this.extendExpansion(kexp, atd, this.initialidx);
            this.Used = this.initialidx;
        end
        
        function extendExpansion(this, kexp, atd, idx)
            % Extends the kernel expansion 'kexp' by the training data
            % sample in 'atd' given by index 'idx'.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % atd: The approximation training data @type data.ApproxTrainData
            % idx: The index of the center to choose @type integer
            kexp.Centers.xi(:,end+1) = atd.xi(:,idx);
            if atd.hasTime
                kexp.Centers.ti(end+1) = atd.ti(idx);
            end
            if atd.hasParams
                kexp.Centers.mui(:,end+1) = atd.mui(:,idx);
            end
            this.Used(end+1) = idx;
        end
        
        function bool = checkStop(this, cnt, rel)
            % Checks the stopping conditions for the adaptive approximation
            % algorithm.
            %
            % Considers maximum expansion size and maximum relative errors
            %
            % See also: MaxExpansionSize MaxRelErr
            vb = KerMor.App.Verbose;
            bool = false;
            if cnt == this.MaxExpansionSize
                if vb > 1
                    fprintf('AAdaptiveBase stopping criteria holds: Max expansion size %d reached.\n',this.MaxExpansionSize);
                end
                bool = true;
            elseif rel < this.MaxRelErr
                if vb > 1
                    fprintf('AAdaptiveBase stopping criteria holds: Relative error %.7e < %.7e\n',rel,this.MaxRelErr);
                end
                bool = true;
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Runs the actual detailed algorithm.
        %
        % Parameters:
        % kexp: The kernel expansion @type kernels.KernelExpansion
        % atd: The approximation training data @type data.ApproxTrainData
        startAdaptiveExtension(this, kexp, atd);
    end
    
    methods(Access=private)
        function [c, idx] = getInitialCenter(this, atd)
            % Computes an initial center from the training data.
            %
            % Depending on the InitialCenter property, the following
            % strategy is pursued:
            % - 'maxfx' The center with the largest (w.r.t. the
            % ABase.ErrorFun norm) associated `f(x_i)` value is chosen.
            % - 'center' The training point closest to the geometrical center
            % if the training data's bounding box is used.
            % - 't0' Tries to find a training sample for time zero and uses
            % that. If None are found, the strategy falls back to 'center'. If
            % multiple are found, the one closest to the geometrical center is
            % chosen.
            %
            % Parameters:
            % atd: The approximation training data @type data.ApproxTrainData
            %
            % Return values:
            % c: The initial center vector `[x; t; \mu]` @type colvec
            % idx: The index of the vector inside the training data @type integer
            %
            % See also: data.ApproxTrainData InitialCenter

            if strcmp(this.InitialCenter,'maxfx')
                [~, idx] = max(this.ErrorFun(atd.fxi));
                c = atd.xi(:,idx);
                if atd.hasTime
                    c = [c; atd.ti(idx)];
                end
                if atd.hasParams
                    c = [c; atd.mui(:,idx)];
                end
            elseif strcmp(this.InitialCenter,'center')
                [c, idx] = atd.getClosestToCenter;
            elseif strcmp(this.InitialCenter,'t0')
                if atd.hasTime
                    idx = find(atd.ti == 0);
                    if ~isempty(idx)
                        c = [atd.xi(:,idx); zeros(1,numel(idx))];
                        if atd.hasParams
                            c = [c; atd.mui(:,idx)];
                        end
                        if numel(idx) > 1
                            [~,i] = min(sum((c-repmat(atd.Center,1,numel(idx))).^2,1));
                            c = c(:,i);
                            idx = idx(i);
                        end
                    end
                end
            else
                error('Unknown strategy: %s',this.InitialCenter);
            end
        end
    end
    
    %% Getter & Setter
    methods
        
        function set.MaxExpansionSize(this, value)
            if ~isposintscalar(value)
                error('Value must be a positive integer.');
            end
            this.MaxExpansionSize = value;
        end
                              
        function set.MaxRelErr(this, value)
            if ~isempty(value) && ~isposrealscalar(value)
                error('The value must be a positive scalar');
            end
            this.MaxRelErr = value;
        end
    end
    
    methods(Static,Access=protected)
        function this = loadobj(this, initfrom)
            if nargin > 1
                this.MaxExpansionSize = initfrom.MaxExpansionSize;
                this.InitialCenter  = initfrom.InitialCenter;
                if isfield(initfrom,'initialidx')
                    this.initialidx = initfrom.initialidx;
                end
                if isfield(initfrom,'Used')
                    this.Used = initfrom.Used;
                end
                this.ExpansionSizes = initfrom.ExpansionSizes;
                this.InitialCenter = initfrom.InitialCenter;
                this.MaxRelErr = initfrom.MaxRelErr;
                this = loadobj@approx.algorithms.ABase(this, initfrom);
            end
        end
    end
end


