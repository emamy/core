classdef ABase < KerMorObject & IParallelizable & ICloneable
% ABase: Base class for any approximation generation algorithms for kernel
% expansions,
%
% @author Daniel Wirtz @date 2011-07-07
%
% @change{0,5,dw,2011-11-02} 
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
% - New default value 'false' for the UsefScaling property; recent experiments suggested 'true'
% might not be a wise default value but an extra source of errors.
%
% @change{0,5,dw,2011-10-14} Improved the parallel computation of
% kernel expansion coefficients.
%
% @change{0,5,dw,2011-09-12} Added initial values that can be passed to
% the CoeffComp algorithms. Now computing coefficients at once if
% MultiTargetComputation of the CoeffComp property is true.
%
% @new{0,5,dw,2011-07-07} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The different kernel expansion configurations to try
        %
        % @propclass{critical} Without this setting this algorithm makes little sense.
        %
        % @type kernels.config.ExpansionConfig @default []
        ExpConfig = [];
        
        % Flag that determines whether the approximation center f values should be scaled to [-1,1]
        % before the approximation is computed.
        %
        % @propclass{optional} This option makes sense when using univariate rotation-invariant
        % kernels as different dimensions might have different scales
        %
        % @default false
        UsefScaling = false;
        
        % The error function to apply on each data vector.
        %
        % Must be a function_handle that takes one matrix argument. The result is computed
        % along the first dimension, i.e. the return value has the same number of columns than
        % the input matrix.
        %
        % @propclass{important} Depending on the approximation goal different error functions
        % are suitable.
        %
        % @type function_handle @default "@Norm.L2"
        %
        % See also: Norm
        ErrorFun = @Norm.L2;
    end
    
    properties(SetAccess=protected)
        % The computation time for the last run in seconds.
        %
        % @type double @default []
        LastCompTime = [];
        
        % For each configuration, contains a row with the maximum errors on the
        % training data.
        % The number of columns depends on the type of algorithm implemented by the subclasses.
        %
        % @default [] @type matrix<double>
        MaxErrors = [];
        
        % For each configuration, contains a row with the maximum relative errors on the
        % training data.
        % The number of columns depends on the type of algorithm implemented by the subclasses.
        %
        % @default [] @type matrix<double>
        MaxRelErrors = [];
    end
    
    methods
        function this = ABase
            this = this@KerMorObject;
            
            this.registerProps('ExpConfig','UsefScaling');
        end
        
        function copy = clone(this, copy)
            copy.ExpConfig = this.ExpConfig;
            copy.UsefScaling = this.UsefScaling;
            copy.ErrorFun = this.ErrorFun;
            copy.MaxErrors = this.MaxErrors;
            copy.MaxRelErrors = this.MaxRelErrors;
            copy.LastCompTime = this.LastCompTime;
        end
        
        
        function computeApproximation(this, kexp, atd)
            time = tic;
            % Scale f-values if wanted
            if this.UsefScaling
                [fm,fM] = atd.fxi.getColBoundingBox;
                s = max(abs(fm),abs(fM));
                s(s==0) = 1;
                oldfxi = atd.fxi.toMemoryMatrix;
                atd.fxi(:,:) = oldfxi ./ repmat(s,1,size(atd.fxi,2));
            end
            
            remove_conf = false;
            if isempty(this.ExpConfig)
                this.ExpConfig = kexp.getDefaultExpansionConfig;
                remove_conf = true;
            else 
                if this.ExpConfig.getNumConfigurations == 0
                    error('Need at least one expansion configuration.');
                end
            end
            
            % Call template method for component wise approximation
            this.templateComputeApproximation(kexp, atd);
            
            % Remove temp config
            if remove_conf
                this.ExpConfig = [];
            end
            
            % Rescale if set
            if this.UsefScaling
                kexp.Ma = diag(s)*kexp.Ma;
                atd.fxi(:,:) = oldfxi;
            end
            
            % Reduce the snapshot array and coeff data to the
            % really used ones! This means if any snapshot x_n is
            % not used in any dimension, it is kicked out at this
            % stage.
            n = size(kexp.Centers.xi,2);
            hlp = sum(abs(kexp.Ma) ~= 0,1);
            usedidx = find(hlp > 0);
            if length(usedidx) < n
                kexp.Ma = kexp.Ma(:,usedidx);
                kexp.Centers.xi = kexp.Centers.xi(:,usedidx);
                if isfield(kexp.Centers,'ti') && ~isempty(kexp.Centers.ti)
                    kexp.Centers.ti = kexp.Centers.ti(:,usedidx);
                end
                if isfield(kexp.Centers,'mui') && ~isempty(kexp.Centers.mui)
                    kexp.Centers.mui = kexp.Centers.mui(:,usedidx);
                end
            end
            
            % Create sparse representation if many elements are zero
            if sum(hlp) / numel(kexp.Ma) < .5
                kexp.Ma = sparse(kexp.Ma);
            end
            this.LastCompTime = toc(time);
        end
    end
    
    methods(Access=protected)
        function [val, idx, errs] = getError(this, kexp, atd)
            % Computes the error according to the chosen error function
            % (see ErrorFun property) with respect to the current kernel
            % expansion and the `f(x_i)` values in the training data.
            %
            % Parameters:
            % kexp: The kernel expansion @type kernels.KernelExpansion
            % atd: The approximation training data @type data.ApproxTrainData
            %
            % Return values:
            % val: The maximum error `\max ||f(x_i,t_i,\mu_i) -
            % \hat{f}(x_i,t_i,\mu_i)||_{\{2,\infty\}}` @type double
            % idx: The index of the maximum error inside the errs vector @type integer
            % errs: A row vector of the errors for each sample @type rowvec
            %
            % See also: ErrorFun
            errs = this.ErrorFun(atd.fxi - kexp.evaluate(atd.xi, atd.ti, atd.mui));
            [val, idx] = max(errs);
        end
    end
    
    methods(Abstract)
        % Plots the errors computed during the last run.
        %
        % Parameters:
        % pm: A PlotManager instance @type PlotManager @default PlotManager
        %
        % Return values:
        % pm: The PlotManager instance, created if none is passed, otherwise the same. @type
        % PlotManager
        pm = plotErrors(this, pm);
    end
    
    methods(Abstract, Access=protected)
        % Performs the actual approximation after scaling.
        %
        % Template method.
        templateComputeApproximation(this, kexp, atd);
    end
    
     methods(Static,Access=protected)
        function this = loadobj(this, initfrom)
            if nargin > 1
                this.ErrorFun = initfrom.ErrorFun;
                this.LastCompTime  = initfrom.LastCompTime;
                this.UsefScaling = initfrom.UsefScaling;
                if isfield(initfrom,'ExpConfig')
                    this.ExpConfig = initfrom.ExpConfig;
                end
                if isfield(initfrom,'MaxErrors') && ~isempty(initfrom.MaxErrors)
                    this.MaxErrors = initfrom.MaxErrors;
                elseif isfield(initfrom,'err')
                    this.MaxErrors = initfrom.err;
                end
                if isfield(initfrom,'MaxRelErrors') && ~isempty(initfrom.MaxRelErrors)
                    this.MaxRelErrors = initfrom.MaxRelErrors;
                elseif isfield(initfrom,'relerr')
                    this.MaxRelErrors = initfrom.relerr;
                end
                this = loadobj@KerMorObject(this, initfrom);
            end
        end
    end
end