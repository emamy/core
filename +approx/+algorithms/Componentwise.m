classdef Componentwise < approx.algorithms.ABase
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
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
%
% @todo think of new structure on how to combine this with the getDists in AAdaptiveBase (or
% extract convenience methods further, general concept of "KernelConfig")
    
    properties(SetObservable)    
        % An instance of a class implementing the approx.algorithms.IKernelCoeffComp
        % interface.
        %
        % This properties class will be used to compute the kernel
        % coefficients for each dimension.
        %
        % @propclass{important} The correct coefficient computation strategy might make the
        % difference.
        %
        % @default general.interpolation.KernelInterpol
        CoeffComp;
        
        % The different coefficient computation algorithm configurations to try
        %
        % @propclass{critical} Without this setting this algorithm makes little sense.
        %
        % @type general.IClassConfig @default []
        CoeffConfig = [];
        
        % Percentage `p` of the training data to use as validation data
        %
        % Admissible values are `p\in]0,\frac{1}{2}]`.
        %
        % @propclass{optional} 
        %
        % @default .2
        % @type double
%         ValidationPercent = .2;
    end
    
    properties(SetAccess=private)
        % The times needed for the coefficients to compute for each exp/coeffcomp configuration.
        %
        % @type matrix<double> @default []
        SingleRuntimes;
    end
            
    methods    
        function this = Componentwise
            this = this@approx.algorithms.ABase;
            
            this.CoeffComp = general.interpolation.KernelInterpol;
            this.CoeffConfig = this.CoeffComp.getDefaultConfig;
            
            % Register default property changed listeners
            this.registerProps('CoeffComp','CoeffConfig');
        end
                        
        function copy = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            copy = approx.algorithms.Componentwise;
            
            copy = clone@approx.algorithms.ABase(this, copy);

            % copy local props
            copy.CoeffComp = this.CoeffComp; % Dont clone the coefficient computation method
            copy.CoeffConfig = this.CoeffConfig;
        end
        
        function pm = plotErrors(this, pm)
            if nargin < 2
                pm = PlotManager(false,1,2);
                pm.LeaveOpen = true;
            end
            
            nc = this.ExpConfig.getNumConfigurations;
            nco = this.CoeffConfig.getNumConfigurations;
            [X,Y] = meshgrid(1:nc,1:nco);
            
            h = pm.nextPlot('abs','Absolute errors','expansion config','coeff comp config');
            ph = tools.LogPlot.logsurf(h,X,Y,this.MaxErrors');
%             set(ph(this.ExpConfig.vBestConfigIndex),'LineWidth',2);
            h = pm.nextPlot('rel','Relative errors','expansion config','coeff comp config');
            ph = tools.LogPlot.logsurf(h,X,Y,this.MaxRelErrors');
%             set(ph(this.ExpConfig.vBestConfigIndex),'LineWidth',2);
            
            if nargin < 2
                pm.done;
            end
        end
    end
    
    methods(Access=protected, Sealed)
        function templateComputeApproximation(this, kexp, atd)
            %
            % @docupdate
             
            % Set centers
            kexp.Centers.xi = atd.xi.toMemoryMatrix;
            kexp.Centers.ti = atd.ti;
            kexp.Centers.mui = atd.mui;
            
            ec = this.ExpConfig;
            nc = ec.getNumConfigurations;
            cc = this.CoeffConfig;
            nco = cc.getNumConfigurations;
            
            % Keep track of maximum errors (matrix-wise for expansion config / coeff config)
            this.MaxErrors = zeros(nc,nco);
            this.MaxRelErrors = zeros(nc,nco);
            this.SingleRuntimes = zeros(nc,nco);
            
            minerr = Inf;
            bestcidx = [];
            bestcoidx = [];
            bestMa = [];
            pi = tools.ProcessIndicator('Trying %d configurations (%d kernel, %d coeffcomp)',...
                nc*nco,false,nc*nco,nc,nco);
            for kcidx = 1:nc
                if KerMor.App.Verbose > 2
                    fprintf('Applying expansion config %s\n',ec.getConfigurationString(kcidx, false));
                end
                ec.applyConfiguration(kcidx, kexp);
                
                % Call coeffcomp preparation method and pass kernel matrix
                this.CoeffComp.init(kexp);

                for coidx = 1:nco
                    if KerMor.App.Verbose > 2
                        fprintf('Applying %s config %s\n',class(this.CoeffComp),cc.getConfigurationString(coidx, false));
                    end
                    cc.applyConfiguration(coidx, this.CoeffComp);
                    
                    % Call protected method
                    time = tic;
                    this.computeCoeffs(kexp, atd.fxi, []);
%                     this.computeCoeffs(kexp, atd.fxi, kexp.Ma);

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

                    if KerMor.App.Verbose > 3
                        fprintf(' ||Ma||:%.5e\n',sum(kexp.Ma_norms));
                    end
                    pi.step;
                end
            end
            
            %% Assign best values
            cc.vBestConfigIndex = bestcoidx;
            cc.applyConfiguration(bestcoidx, this.CoeffComp);
            ec.setBestConfig(bestcidx, kexp);
            kexp.Ma = bestMa;
            
            if KerMor.App.Verbose > 1
                figure;
                plot(this.MaxErrors,'r');
            end
            
            pi.stop;
        end
    end
    
    methods(Access=private)
    
        function computeCoeffs(this, kexp, fxi, initialalpha)
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
            if nargin < 4 || isempty(initialalpha)
                initialalpha = double.empty(size(fxi,1),0);
            end
            if this.ComputeParallel
                if KerMor.App.Verbose > 3
                    fprintf('Starting parallel component-wise coefficient computation\n');
                end
                this.computeCoeffsParallel(kexp, fxi, initialalpha);
            else
                if KerMor.App.Verbose > 3
                    fprintf('Starting component-wise coefficient computation\n');
                end
                this.computeCoeffsSerial(kexp, fxi, initialalpha);
            end
        end
    
        
        function computeCoeffsSerial(this, kexp, fxi, initialalpha)
            % Computes the coefficients using the CoeffComp instance
            % serially.
            %
            % @throws KerMor:coeffcomp_failed Forwarded from CoeffComp
            %
            % @todo remove waitbar and connect to verbose/messaging system
            %% Non-parallel execution
            fdims = size(fxi,1);
            n = size(kexp.Centers.xi,2);
            kexp.Ma = zeros(fdims, n);
            if this.CoeffComp.MultiTargetComputation
                if KerMor.App.Verbose > 3
                    fprintf('Computing approximation for all %d dimensions...\n',fdims);
                end
                % Call template method
                [ai, svidx] = this.CoeffComp.computeKernelCoefficients(fxi,initialalpha);
                kexp.Ma(:,svidx) = ai;
            else
                for fdim = 1:fdims
                    if KerMor.App.Verbose > 3 && fdims > 1
                        fprintf('Computing approximation for dimension %d/%d ... %2.0f %%\n',fdim,fdims,(fdim/fdims)*100);
                    end
                    % Call template method
                    [ai, svidx] = this.CoeffComp.computeKernelCoefficients(fxi(fdim,:),initialalpha(fdim,:)); 
                    kexp.Ma(fdim,svidx) = ai;
                end
            end
        end
        
        function computeCoeffsParallel(this, kexp, fxi, initialalpha)
            %% Parallel execution
            n = size(kexp.Centers.xi,2);
            fdims = size(fxi,1);
            fprintf('Starting parallel component-wise approximation computation of %d dimensions on %d workers...\n',fdims,matlabpool('size'));
            Ma = zeros(fdims, n);
            parfor fdim = 1:fdims
                % Call template method
                [ai, sv] = this.CoeffComp.computeKernelCoefficients(...
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
        function set.CoeffComp(this, value)
            if ~isa(value,'approx.algorithms.IKernelCoeffComp')
                error('Property value must implement the approx.algorithms.IKernelCoeffComp interface');
            end
            this.CoeffComp = value;
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
            end
        end
    end
end


