classdef BaseCompWiseKernelApprox < approx.BaseApprox & ...
        dscomponents.CompwiseKernelCoreFun & IParallelizable %& approx.IAutoConfig
    %Base class for component-wise kernel approximations
    %
    % For each dimension `k` there is a representation
    % ``f_k(x) = \sum\limits_{i=1}^N \alpha_{k,i}\Phi(x,x_i)``
    % for the approximation. The property cData contains in row `k` all
    % indices `\alpha_k` used in the `k`-th dimension.
    %
    % The property Centers (inherited from
    % dscomponents.CompwiseKernelCoreFun) contains all state variable
    % Centers that are relevant for the evaluation of the approximated
    % function. (No matter how many originally have been used for the
    % approximation computation)
    %
    % @change{0,4,dw,2011-05-03} Removed off-property dependencies in computeCoeffs routines as
    % expansion offsets are no longer used.
    %
    % @change{0,3,dw,2011-04-26} Implementing the approximateCoreFun in this class so that rescaling
    % of the function f values `f_{x_i}` can be used. Subclasses now implement the template method
    % approx.BaseCompWiseKernelApprox.computeCompwiseApprox.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @change{0,3,dw,2011-03-31} 
    % - Added a new property approx.BaseCompWiseKernelApprox.CoeffComp
    % which allows to use the strategy pattern for computing kernel
    % expansion coefficients.
    % - Moved some old methods into the new class
    % approx.DefaultCompWiseKernelApprox to retain the old functionality
    % which selects a subset of the projection training data by linspace
    % selection as approximation data.
    % - Changed method names for computeSerial/Parallel to
    % computeCoeffsSerial/computeCoeffsParallel and visibility to protected
    % as they can be used by any subclass.
    %
    % See also: BaseKernelApprox

    properties(SetObservable)
        % An instance of a class implementing the approx.IKernelCoeffComp
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
        
        function this = BaseCompWiseKernelApprox
            this = this@approx.BaseApprox;
            this = this@dscomponents.CompwiseKernelCoreFun;
            
            this.CoeffComp = general.interpolation.KernelInterpol;
            
            this.registerProps('CoeffComp','UsefScaling');
        end
        
        function approximateCoreFun(this, model)
            % 
            % @todo find out when sparse representation is more
            % efficient!
            
            fxi = model.Data.ApproxfValues;
            % Scale f-values if wanted
            if this.UsefScaling
                [fm,fM] = general.Utils.getBoundingBox(fxi);
                s = max(abs(fm),abs(fM));
                s(s==0) = 1;
                fxi = fxi ./ repmat(s,1,size(fxi,2));
            end
            
            % Call template method for component wise approximation
            this.computeCompwiseApprox(model,fxi);
            
            % Rescale if set
            if this.UsefScaling
                this.Ma = diag(s)*this.Ma;
            end
            
            % Reduce the snapshot array and coeff data to the
            % really used ones! This means if any snapshot x_n is
            % not used in any dimension, it is kicked out at this
            % stage.
            n = size(this.Centers.xi,2);
            hlp = sum(this.Ma ~= 0,1);
            usedidx = find(hlp > 0);
            if length(usedidx) < n
                this.Ma = this.Ma(:,usedidx);
                this.Centers.xi = this.Centers.xi(:,usedidx);
                if ~isempty(this.Centers.ti)
                    this.Centers.ti = this.Centers.ti(:,usedidx);
                end
                if ~isempty(this.Centers.mui)
                    this.Centers.mui = this.Centers.mui(:,usedidx);
                end
            end
            
            % Create sparse representation if many elements are zero
            if sum(hlp) / numel(this.Ma) < .5
                this.Ma = sparse(this.Ma);
            end
        end
        
        function target = clone(this, target)
            % Clones the instance.
            %
            % Note:
            % Since we use multiple inheritance we have to call the clone
            % method from both above superclasses. In this case this leads
            % to a double execution of the clone method from
            % dscomponents.ACoreFcn, but this is rather done than ommitted
            % and causing trouble later on.
            if nargin == 1 || ~(isa(target,'approx.BaseApprox') && isa(target,'dscomponents.CompwiseKernelCoreFun'))
                error('Invalid clone call. As this class is abstract, a subclass of both BaseApprox and CompwiseKernelCoreFun has to be passed as second argument.');
            end
                       
            target = clone@dscomponents.CompwiseKernelCoreFun(this, target);
            target = clone@approx.BaseApprox(this, target);
            
            % copy local props
            copy.CoeffComp = this.CoeffComp;
            copy.UsefScaling = this.UsefScaling;
        end
        
        function set.CoeffComp(this, value)
            if ~isa(value,'approx.IKernelCoeffComp')
                error('Property value must implement the approx.IKernelCoeffComp interface');
            end
            this.CoeffComp = value;
        end
    end
    
    methods(Access=protected)
        function computeCoeffs(this, fxi)
            % Computes the coefficients for all components.
            %
            % Convenience method that any subclasses may (must?) use in
            % order to compute the coefficients. Depending on the parallel
            % flag this is done parallel or not (parfor)
            %
            % Please take care that the CoeffComp.init method was called
            % before executing this function.
%             if KerMor.App.Verbose > 1
%                 fprintf('Initializing component-wise coefficient computation\n');
%             end
            if this.ComputeParallel
                this.computeCoeffsParallel(fxi);
            else
                this.computeCoeffsSerial(fxi);
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Inner approximation computation
        %
        % Instead of accessing model.Data.ApproxfValues those value are pre-processed if scaling is
        % used and thus passed as argument to the inner approximation computation.
        %
        % Parameters:
        % model: The full model
        % fxi: The full core function f evaluations at the ApproxTrainData centers. Either scaled or
        % not scaled.
        computeCompwiseApprox(this, model, fxi);
    end
    
    methods(Access=private)
        
        function computeCoeffsSerial(this, fxi)
            % Computes the coefficients using the CoeffComp instance
            % serially.
            %
            % @todo remove waitbar and connect to verbose/messaging system
            %% Non-parallel execution
            n = size(this.Centers.xi,2);
            fdims = size(fxi,1);
            this.Ma = zeros(fdims, n);
            for fdim = 1:fdims
                if KerMor.App.Verbose > 3
                    fprintf('Computing approximation for dimension %d/%d ... %2.0f %%\n',fdim,fdims,(fdim/fdims)*100);
                end
                % Call template method
                [ai, svidx] = this.CoeffComp.computeKernelCoefficients(fxi(fdim,:));
                if ~isempty(svidx)
                    this.Ma(fdim,svidx) = ai;
                else
                    this.Ma(fdim,:) = ai;
                end
            end
        end
        
        function computeCoeffsParallel(this, fxi)
            %% Parallel execution
            n = size(this.Centers.xi,2);
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
            
            this.Ma = zeros(fdims, n);
            for idx = 1:fdims
                if ~isempty(parSV{idx})
                    this.Ma(idx,parSV{idx}) = parAI{idx};
                else
                    this.Ma(idx,:) = parAI{idx};
                end
            end
        end
    end
    
end


