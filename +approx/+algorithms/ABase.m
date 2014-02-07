classdef ABase < KerMorObject & ICloneable & IReductionSummaryPlotProvider
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        % @type function_handle @default Norm.L2
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
        
        % For each effective configuration, the stop flags are stored here.
        %
        % @default [] @type matrix<integer>
        %
        % See also: StopFlag
        StopFlags = [];
        
        % Index of the best expansion config determined by the algorithm
        %
        % @type integer
        BestExpConfig;
    end
    
    properties(GetAccess=protected,SetAccess=private)
        % If UsefScaling is set to true, this matrix contains the scaling matrix
        % which has to be used in sub-algorithms in order to compute the correct approximation
        % error on training data.
        %
        % This is automatically done upon calls to ABase.getError.
        %
        % See also: UsefScaling
        ScalingG;
    end
    
    methods
        function this = ABase
            this = this@KerMorObject;
            
            this.registerProps('ExpConfig','UsefScaling');
        end
        
        function copy = clone(this, copy)
            if ~isempty(this.ExpConfig)
                copy.ExpConfig = this.ExpConfig.clone;
            end
            copy.UsefScaling = this.UsefScaling;
            copy.ErrorFun = this.ErrorFun;
            copy.MaxErrors = this.MaxErrors;
            copy.MaxRelErrors = this.MaxRelErrors;
            copy.LastCompTime = this.LastCompTime;
            copy.StopFlags = this.StopFlags;
            copy.ScalingG = this.ScalingG;
            copy.BestExpConfig = this.BestExpConfig;
        end
        
        
        function kexp = computeApproximation(this, atd)
            time = tic;
            
            if isempty(this.ExpConfig)
                error('No ExpConfig set. Aborting.');
            elseif this.ExpConfig.getNumConfigurations == 0
                error('Need at least one expansion configuration.');
            end
            
            % Scale f-values if wanted
            if this.UsefScaling
                [fm,fM] = atd.fxi.getColBoundingBox;
                s = max(abs(fm),abs(fM));
                s(s==0) = 1;
                oldfxi = atd.fxi.toMemoryMatrix;
                atd.fxi(:,:) = oldfxi ./ repmat(s,1,size(atd.fxi,2));
                this.ScalingG = diag(s);
            else
                this.ScalingG = 1;
            end
            
            % Call template method for component wise approximation
            kexp = this.templateComputeApproximation(atd);
            
            % Rescale if set
            if this.UsefScaling
                kexp.Ma = this.ScalingG*kexp.Ma;
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
        
        function [str, rangetab] = getApproximationSummary(this)
            % Setup
            str = [object2str(this,1) char(13)];
            
            ec = this.ExpConfig;
            % Computation time
            str = [str sprintf('Total computation time for %d configurations: %gs (%gmin,%gh)\n',...
                ec.getNumConfigurations,this.LastCompTime,...
                this.LastCompTime/60,this.LastCompTime/3600)];
            % Configuration
            str = [str sprintf('Best expansion configuration at index %d:\n%s\n',this.BestExpConfig,...
                ec.getConfigurationString(this.BestExpConfig))];
            rangetab = ec.getValueRanges;
            if nargout < 2
                str = [str sprintf('Expansion configuration ranges:\n%s\n',...
                    rangetab.toString)];
                if nargout < 1
                    disp(str);
                end
            end
        end
        
        function plotSummary(this, pm, context)
            if nargin < 3
                context = 'Summary';
                if nargin < 2
                    pm = PlotManager(false,1,2);
                    pm.LeaveOpen = true;
                end
            end
            
            str = sprintf('%s: Max absolute errors',context);
            h = pm.nextPlot('maxerrors',str,...
                'expansion size','error');
            ph = semilogy(h,1:size(this.MaxErrors,2),this.MaxErrors');
            if ~isempty(this.ExpConfig)
                set(ph(this.BestExpConfig),'LineWidth',2);
            end
            str = sprintf('%s: Max relative errors',context);
            h = pm.nextPlot('maxrelerrors',str,...
                'expansion size','error');
            ph = semilogy(h,1:size(this.MaxRelErrors,2),this.MaxRelErrors');
            if ~isempty(this.ExpConfig)
                set(ph(this.BestExpConfig),'LineWidth',2);
            end
            
            if nargin < 2
                pm.done;
            end
        end
        
        function nc = getTotalNumConfigurations(this)
            nc = 1;
            if ~isempty(this.ExpConfig)
                nc = this.ExpConfig.StateConfig.getNumConfigurations;
            end
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
            e = atd.fxi - kexp.evaluate(atd.xi, atd.ti, atd.mui);
            errs = this.ErrorFun(this.ScalingG*e);
            [val, idx] = max(errs);
        end
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
                if ~isfield(initfrom,'StopFlags') || isempty(initfrom.StopFlags)
                    this.StopFlags = zeros(size(this.MaxErrors));
                end
                this = loadobj@KerMorObject(this, initfrom);
            else
                if isempty(this.StopFlags)
                    this.StopFlags = zeros(size(this.MaxErrors));
                end
                this = loadobj@KerMorObject(this);
            end
        end
    end
end