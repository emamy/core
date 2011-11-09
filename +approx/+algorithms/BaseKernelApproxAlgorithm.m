classdef BaseKernelApproxAlgorithm < KerMorObject & IParallelizable & ICloneable
% BaseKernelApproxAlgorithm: Base class for any approximation generation algorithms for kernel
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
        
        % Flag that determines whether the approximation center f values should be scaled to [-1,1]
        % before the approximation is computed.
        %
        % @propclass{optional} This option makes sense when using univariate rotation-invariant
        % kernels as different dimensions might have different scales
        %
        % @default false
        UsefScaling = false;
    end
    
    methods
        function this = BaseKernelApproxAlgorithm
            this = this@KerMorObject;
            
            this.CoeffComp = general.interpolation.KernelInterpol;
            
            this.registerProps('CoeffComp','UsefScaling');
        end
        
        function copy = clone(this, copy)
            copy.CoeffComp = this.CoeffComp; % Dont clone the coefficient computation method
            copy.UsefScaling = this.UsefScaling;
        end
        
        
        function computeApproximation(this, kexp, atd)
            
            % Scale f-values if wanted
            if this.UsefScaling
                [fm,fM] = general.Utils.getBoundingBox(atd.fxi);
                s = max(abs(fm),abs(fM));
                s(s==0) = 1;
                oldfxi = atd.fxi;
                atd.fxi = atd.fxi ./ repmat(s,1,size(atd.fxi,2));
            end
            
            % Call template method for component wise approximation
            this.detailedComputeApproximation(kexp, atd);
            
            % Rescale if set
            if this.UsefScaling
                kexp.Ma = diag(s)*kexp.Ma;
                atd.fxi = oldfxi;
            end
            
            % Reduce the snapshot array and coeff data to the
            % really used ones! This means if any snapshot x_n is
            % not used in any dimension, it is kicked out at this
            % stage.
            n = size(kexp.Centers.xi,2);
            hlp = sum(kexp.Ma ~= 0,1);
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
    
    methods(Access=protected)
        function computeCoeffs(this, kexp, fxi, initialalpha)
            % Computes the coefficients for all components.
            %
            % Convenience method that any subclasses may (must?) use in
            % order to compute the coefficients. Depending on the parallel
            % flag this is done parallel or not (parfor)
            %
            % Please take care that the CoeffComp.init method was called
            % before executing this function.
            if isempty(initialalpha)
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
    end
    
    methods(Access=private)
        
        function computeCoeffsSerial(this, kexp, fxi, initialalpha)
            % Computes the coefficients using the CoeffComp instance
            % serially.
            %
            % @throws KerMor:coeffcomp_failed Forwarded from CoeffComp
            %
            % @todo remove waitbar and connect to verbose/messaging system
            %% Non-parallel execution
            fdims = size(fxi,1);
            oldma = kexp.Ma;
            try
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
                        if KerMor.App.Verbose > 3
                            fprintf('Computing approximation for dimension %d/%d ... %2.0f %%\n',fdim,fdims,(fdim/fdims)*100);
                        end
                        % Call template method
                        [ai, svidx] = this.CoeffComp.computeKernelCoefficients(fxi(fdim,:),initialalpha(fdim,:)); 
                        kexp.Ma(fdim,svidx) = ai;
                    end
                end
            catch ME
                if strcmp(ME.identifier,'KerMor:coeffcomp_failed')
                    % Restore old coefficients
                    kexp.Ma = oldma;
                end
                rethrow(ME);
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
    
    methods(Abstract, Access=protected)
        % Performs the actual approximation after scaling.
        %
        % Template method.
        detailedComputeApproximation(this, kexp, atd);
    end
    
end