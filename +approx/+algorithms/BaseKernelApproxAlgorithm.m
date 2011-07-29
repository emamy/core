classdef BaseKernelApproxAlgorithm < KerMorObject & IParallelizable
    % BaseKernelApproxAlgorithm: Base class for any approximation generation algorithms for kernel
    % expansions,
    %
    % @author Daniel Wirtz @date 2011-07-07
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
        % @default true
        UsefScaling = true;
    end
    
    methods
        function this = BaseKernelApproxAlgorithm
            this = this@KerMorObject;
            
            this.CoeffComp = general.interpolation.KernelInterpol;
            
            this.registerProps('CoeffComp','UsefScaling');
        end
        
        
        function computeApproximation(this, kexp, xi, ti, mui, fxi)
            
            % Scale f-values if wanted
            if this.UsefScaling
                [fm,fM] = general.Utils.getBoundingBox(fxi);
                s = max(abs(fm),abs(fM));
                s(s==0) = 1;
                fxi = fxi ./ repmat(s,1,size(fxi,2));
            end
            
            % Call template method for component wise approximation
            this.detailedComputeApproximation(kexp, xi, ti, mui, fxi);
            
            % Rescale if set
            if this.UsefScaling
                kexp.Ma = diag(s)*kexp.Ma;
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
                if ~isempty(kexp.Centers.ti)
                    kexp.Centers.ti = kexp.Centers.ti(:,usedidx);
                end
                if ~isempty(kexp.Centers.mui)
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
        function computeCoeffs(this, kexp, fxi)
            % Computes the coefficients for all components.
            %
            % Convenience method that any subclasses may (must?) use in
            % order to compute the coefficients. Depending on the parallel
            % flag this is done parallel or not (parfor)
            %
            % Please take care that the CoeffComp.init method was called
            % before executing this function.
            if this.ComputeParallel
                if KerMor.App.Verbose > 2
                    fprintf('Starting parallel component-wise coefficient computation\n');
                end
                this.computeCoeffsParallel(kexp, fxi);
            else
                if KerMor.App.Verbose > 2
                    fprintf('Starting component-wise coefficient computation\n');
                end
                this.computeCoeffsSerial(kexp, fxi);
            end
        end
    end
    
    methods(Access=private)
        
        function computeCoeffsSerial(this, kexp, fxi)
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
                for fdim = 1:fdims
                    if KerMor.App.Verbose > 3
                        fprintf('Computing approximation for dimension %d/%d ... %2.0f %%\n',fdim,fdims,(fdim/fdims)*100);
                    end
                    % Call template method
                    [ai, svidx] = this.CoeffComp.computeKernelCoefficients(fxi(fdim,:)); 
                    if ~isempty(svidx)
                        kexp.Ma(fdim,svidx) = ai;
                    else
                        kexp.Ma(fdim,:) = ai;
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
        
        function computeCoeffsParallel(this, kexp, fxi)
            %% Parallel execution
            n = size(kexp.Centers.xi,2);
            fdims = size(fxi,1);
            fprintf('Starting parallel component-wise approximation computation of %d dimensions on %d workers...\n',fdims,matlabpool('size'));
            parAI = cell(fdims, n);
            parSV = cell(1,fdims);
            % Create handle to speedup communications (otherwise the
            % whole object will be copied to each worker)
            cfun = @this.CoeffComp.computeKernelCoefficients;
            parfor fdim = 1:fdims
                %waitbar(fdim/fdims+10,wh,sprintf('Computing
                %approximation for dimension %d/%d ... %2.0f %%',fdim,fdims,(fdim/fdims)*100));
                % Call template method
                [parAI{fdim}, parSV{fdim}] = cfun(fxi(fdim,:));
            end
            
            kexp.Ma = zeros(fdims, n);
            for idx = 1:fdims
                if ~isempty(parSV{idx})
                    kexp.Ma(idx,parSV{idx}) = parAI{idx};
                else
                    kexp.Ma(idx,:) = parAI{idx};
                end
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Performs the actual approximation after scaling.
        %
        % Template method.
        detailedComputeApproximation(this, kexp, xi, ti, mui, fxi);
    end
    
end