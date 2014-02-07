classdef Componentwise < approx.algorithms.ABase & IParallelizable
% Componentwise: Component-wise kernel approximation with fixed center set.
%
% @author Daniel Wirtz @date 2011-06-01
%
% @change{0,7,dw,2011-11-22} Moved Componentwise.guessGammas to
% kernels.config.ExpansionConfig
%
% @change{0,5,dw,2011-11-02} 
% - New interface for approximation computation: Passing an data.ApproxTrainData instance now
% instead of 'xi,ti,mui' parameters.
% - Re-enabled the clone method
% - Using the new data.ApproxTrainData for guessGammas
%
% @new{0,5,dw,2011-10-14}
% - Added the method Componentwise.guessGammas as a helper method to determine suitable
% Gaussian kernel configurations.
% - This algorithm now also works with kernels.ParamTimeKernelExpansion's
% 
% @new{0,5,dw,2011-07-07} Moved the old approx.Componentwise class to this class.
%
% @new{0,4,dw,2011-06-01} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
%
% @todo think of new structure on how to combine this with the getDists in AAdaptiveBase (or
% extract convenience methods further, general concept of "KernelConfig")
    
    properties(SetObservable)    
        % The different coefficient computation algorithm configurations to
        % try. The IClassConfig.Prototype is used as actual algorithm, and
        % needs to implement the IKernelCoeffComp interface.
        %
        % @propclass{critical} Without this setting this algorithm makes
        % little sense.
        %
        % @type IClassConfig @default general.interpolation.InterpolConfig
        CoeffConfig = [];
        
        % The index of the best coefficient computation configuration.
        %
        % @propclass{data} This property is the result if algorithm
        % execution.
        %
        % @type integer @default []
        %
        % See also: CoeffConfig
        BestCoeffConfig;
    end
    
    properties(SetAccess=private)
        % The times needed for the coefficients to compute for each exp/coeffcomp configuration.
        %
        % @type matrix<double> @default []
        SingleRuntimes;
    end
    
    properties(Access=private, Transient)
        % A temporarily stored IKernelCoeffComp instance
        ccomp;
    end
            
    methods    
        function this = Componentwise
            this = this@approx.algorithms.ABase;
            
            this.CoeffConfig = general.interpolation.InterpolConfig;
            
            % Register default property changed listeners
            this.registerProps('CoeffConfig');
        end
                        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.Componentwise;
            copy = clone@approx.algorithms.ABase(this, copy);

            % copy local props
            copy.BestCoeffConfig = this.BestCoeffConfig;
            if ~isempty(this.CoeffConfig)
                copy.CoeffConfig = this.CoeffConfig.clone;
            end
            copy.SingleRuntimes = this.SingleRuntimes;
        end
        
        function pm = plotSummary(this, pm)
            % Overrides the approx.algorithms.ABase function
            if nargin < 2
                pm = PlotManager(false,1,2);
                pm.LeaveOpen = true;
            end
            
            nc = this.ExpConfig.getNumConfigurations;
            nco = this.CoeffConfig.getNumConfigurations;
            [X,Y] = meshgrid(1:nc,1:nco);
            
            h = pm.nextPlot('abs','Absolute errors','expansion config','coeff comp config');
            LogPlot.logsurf(h,X,Y,this.MaxErrors');
            h = pm.nextPlot('rel','Relative errors','expansion config','coeff comp config');
            LogPlot.logsurf(h,X,Y,this.MaxRelErrors');
            if nargin < 2
                pm.done;
            end
        end
        
        function [str, rangetab] = getApproximationSummary(this)
            [str, rangetab] = getApproximationSummary@approx.algorithms.ABase(this);
            cc = this.CoeffConfig;
            str = [str sprintf('Best coefficient comp configuration at index %d:\n%s\n',this.BestCoeffConfig,...
                cc.getConfigurationString(this.BestCoeffConfig))];
            rangetab = rangetab.append(cc.getValueRanges);
            if nargout < 2
                str = [str sprintf('Configuration ranges:\n%s\n',...
                    rangetab.toString)];
                if nargout < 1
                    disp(str);
                end
            end
        end
        
        function nc = getTotalNumConfigurations(this)
            nc = getTotalNumConfigurations@approx.algorithms.ABase(this);
            if ~isempty(this.CoeffConfig)
                nc = nc * this.CoeffConfig.getNumConfigurations;
            end
        end
    end
    
    methods(Access=protected, Sealed)
        function kexp = templateComputeApproximation(this, atd)
            ec = this.ExpConfig;
            
            kexp = ec.Prototype;
            kexp.Centers.xi = atd.xi.toMemoryMatrix;
            kexp.Centers.ti = atd.ti;
            kexp.Centers.mui = atd.mui;
                
            nc = ec.getNumConfigurations;
            cc = this.CoeffConfig;
            nco = cc.getNumConfigurations;
            
            if nco == 0
                error('No coefficient computation configurations set. See CoeffConfig');
            end
            
            % Keep track of maximum errors (matrix-wise for expansion config / coeff config)
            this.MaxErrors = zeros(nc,nco);
            this.MaxRelErrors = zeros(nc,nco);
            this.StopFlags = zeros(nc,nco);
            this.SingleRuntimes = zeros(nc,nco);
            
            minerr = Inf;
            bestcidx = [];
            bestcoidx = [];
            bestMa = [];
            pi = ProcessIndicator('Trying %d configurations (%d kernel, %d coeffcomp)',...
                nc*nco,false,nc*nco,nc,nco);
            for kcidx = 1:nc
                if KerMor.App.Verbose > 2
                    fprintf('Applying expansion config %s\n',ec.getConfigurationString(kcidx, false));
                end
                kexp = ec.configureInstance(kcidx);
                
                % Initialize the prototype for coeffconfigs with the
                % current kernel expansion
                cc.Prototype.init(kexp);

                for coidx = 1:nco
                    
                    this.ccomp = cc.configureInstance(coidx);
                    
                    if KerMor.App.Verbose > 2
                        fprintf('Applying %s config %s\n',class(this.ccomp),cc.getConfigurationString(coidx, false));
                    end
                    
                    % Call protected method
                    time = tic;
                    sf = this.computeCoeffs(kexp, atd.fxi, []);
%                     sf = this.computeCoeffs(kexp, atd.fxi, kexp.Ma);
                    this.StopFlags(kcidx,coidx) = sf;

                    % Determine maximum error
                    [val, maxidx] = this.getError(kexp, atd);
                    rel = val / (norm(atd.fxi(:,maxidx))+eps);
                    this.MaxErrors(kcidx,coidx) = val;
                    this.MaxRelErrors(kcidx,coidx) = rel;
                    this.SingleRuntimes(kcidx,coidx) = toc(time);
                    
                    if val < minerr
                        minerr = val;
                        bestMa = kexp.Ma;
                        bestcidx = kcidx;
                        bestcoidx = coidx;
                        if KerMor.App.Verbose > 3
                            fprintf(' b: %.5e (rel %.5e)',val,rel);
                        end
                    else
                        if KerMor.App.Verbose > 3
                            fprintf(' w: %.5e  (rel %.5e)',val,rel);
                        end
                    end
                    pi.step;
                end
            end
            
            %% Assign best values
            this.BestCoeffConfig = bestcoidx;
            this.BestExpConfig = bestcidx;
            kexp = ec.configureInstance(bestcidx);
            kexp.Ma = bestMa;
            
            if KerMor.App.Verbose > 1
                figure;
                plot(this.MaxErrors,'r');
            end
            
            pi.stop;
        end
    end
    
    methods(Access=private)
    
        function sf = computeCoeffs(this, kexp, fxi, initialalpha)
            % Computes the coefficients for all components.
            %
            % Convenience method that any subclasses may (must?) use in
            % order to compute the coefficients. Depending on the parallel
            % flag this is done parallel or not (parfor)
            %
            % Please take care that the CoeffComp.init method was called
            % before executing this function.
            %
            % Parameters:
            % kexp: The kernel expansion
            % fxi: The `f(x_i)`values at the expansion centers
            % initialalpha: Initial `\alpha_i` value to use as
            % initialization (if applicable for the algorithm)
            %
            % Return values:
            % sf: The stop flag. @type integer
            %
            % See also: StopFlag
            if nargin < 4 || isempty(initialalpha)
                initialalpha = double.empty(size(fxi,1),0);
            end
            if this.ComputeParallel
                if KerMor.App.Verbose > 3
                    fprintf('Starting parallel component-wise coefficient computation\n');
                end
                sf = this.computeCoeffsParallel(kexp, fxi, initialalpha);
            else
                if KerMor.App.Verbose > 3
                    fprintf('Starting component-wise coefficient computation\n');
                end
                sf = this.computeCoeffsSerial(kexp, fxi, initialalpha);
            end
        end
        
        function sf = computeCoeffsSerial(this, kexp, fxi, initialalpha)
            % Computes the coefficients using the CoeffComp instance
            % serially.
            %
            %
            % @todo remove waitbar and connect to verbose/messaging system
            %
            % Return values:
            % sf: The stop flag. @type integer
            %
            % See also: StopFlag
            
            fdims = size(fxi,1);
            n = size(kexp.Centers.xi,2);
            kexp.Ma = zeros(fdims, n);
            if this.ccomp.MultiTargetComputation
                if KerMor.App.Verbose > 3
                    fprintf('Computing approximation for all %d dimensions...\n',fdims);
                end
                % Call template method
                [ai, svidx, sf] = this.ccomp.computeKernelCoefficients(fxi,initialalpha);
                kexp.Ma(:,svidx) = ai;
            else
                for fdim = 1:fdims
                    if KerMor.App.Verbose > 3 && fdims > 1
                        fprintf('Computing approximation for dimension %d/%d ... %2.0f %%\n',fdim,fdims,(fdim/fdims)*100);
                    end
                    % Call template method
                    [ai, svidx, sf] = this.ccomp.computeKernelCoefficients(fxi(fdim,:),initialalpha(fdim,:)); 
                    kexp.Ma(fdim,svidx) = ai;
                end
            end
        end
        
        function sf = computeCoeffsParallel(this, kexp, fxi, initialalpha)
            % Parallel execution
            %
            % Return values:
            % sf: The stop flag. @type integer
            %
            % See also: StopFlag
            n = size(kexp.Centers.xi,2);
            fdims = size(fxi,1);
            fprintf('Starting parallel component-wise approximation computation of %d dimensions on %d workers...\n',fdims,matlabpool('size'));
            Ma = zeros(fdims, n);
            parfor fdim = 1:fdims
                % Call template method
                [ai, sv, sf] = this.ccomp.computeKernelCoefficients(...
                    fxi(fdim,:),initialalpha(fdim,:));%#ok
                Ai = zeros(1,n);
                Ai(sv) = ai;
                Ma(fdim,:) = Ai;
            end
            kexp.Ma = Ma;
        end
    end
    
    %% Getter & Setter
    methods
        function set.CoeffConfig(this, value)
            if ~isempty(value) && ~isa(value,'IClassConfig')
                error('Property value must implement the IClassConfig interface');
            end
            this.CoeffConfig = value;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this,'approx.algorithms.Componentwise')
                newinst = approx.algorithms.Componentwise;
                if isfield(this, 'SingleRuntimes')
                    newinst.SingleRuntimes = this.SingleRuntimes;
                end
                this = loadobj@approx.algorithms.ABase(newinst, this);
            else
                this = loadobj@approx.algorithms.ABase(this);
            end
        end
    end
    
    methods(Static)
        function res = test_Componentwise
            m = models.pcd.PCDModel(1);
            m.dt = .1;
            m.T = 1;
            m.EnableTrajectoryCaching = false;
            
            % Define prototype expansion
            ec = kernels.config.ExpansionConfig;
            ec.StateConfig = kernels.config.GaussConfig('G',.1:.3);
            
            a = approx.algorithms.Componentwise;
            a.ExpConfig = ec;
            a.CoeffConfig = general.interpolation.InterpolConfig;
            
            ap = approx.KernelApprox;
            ap.Algorithm = a;
            
            m.Approx = ap;
            m.System.Params(1).Desired = 2;
            m.SpaceReducer = spacereduction.PODGreedy;
            m.SpaceReducer.Eps = 1e-2;
            m.offlineGenerations;
            
            mu = m.getRandomParam;
            r = m.buildReducedModel;
            r.simulate(mu);
            
            a.CoeffConfig = general.regression.EpsSVRConfig([.1 .15; 1 2]);
            svr = general.regression.ScalarEpsSVR_SMO;
            svr.MaxCount = 20;
            a.CoeffConfig.Prototype = svr;
            a.UsefScaling = true;
            m.off5_computeApproximation;
            
            r = m.buildReducedModel;
            r.simulate(mu);

            res = true;
        end
    end
end


