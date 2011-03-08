classdef BaseCompWiseKernelApprox < approx.BaseApprox & ...
        dscomponents.CompwiseKernelCoreFun & IParallelizable
    %Base class for component-wise kernel approximations
    %
    % For each dimension `k` there is a representation
    % ``f_k(x) = \sum\limits_{i=1}^N \alpha_{k,i}\Phi(x,x_i) + b_k``
    % for the approximation. The property cData contains in row `k` all
    % indices `\alpha_k` used in the `k`-th dimension. off contains all
    % `b_k`.
    %
    % The property snData contains all state variable snData that
    % are relevant for the evaluation of the approximated function. (No
    % matter how many originally have been used for the approximation
    % computation!)
    %
    % See also: BaseKernelApprox

    properties
       % The number of projection training data snapshots used to compile
        % the approximation training data set. So far, the default strategy
        % implemented in this class simply uses linspace to select a subset
        % of the specified size.
        %
        % Default: 120
        ApproxExpansionSize = 120;
    end
    
    methods
        
        function guessKernelConfig(this,model)
            zero = 0.01;
            trange = 5; srange = 5;
            %prange = 2;
            
            atd = model.Data.ApproxTrainData;
            xd = sqrt(sum(atd(4:end,:).^2));
            v = unique(round(atd(1,:)));
            maxdiff = zeros(1,length(v));
            for muidx = 1:length(v)
                sel = atd(1,:) == v(muidx);
                tmp = xd(sel);
                maxdiff(muidx) = max(abs(tmp(1:end-1)-tmp(2:end)));
            end
            d = max(maxdiff);
            sgamma = -((srange*d)^2)/log(zero);
            
            %             params = model.Data.getParams(atd(1,:));
            %             mud = sqrt(sum(.^2));
            %
            %             pgamma = -((prange*d)^2)/log(zero);
            tgamma = -(trange^2*model.dt)/log(zero);
            
            %fprintf('Setting kernel gammas: Time:%e, System:%e, Params:%e\n',tgamma,sgamma,pgamma);
            fprintf('Setting kernel gammas: Time:%10.10f, System:%10.10f\n',tgamma,sgamma);
            this.TimeKernel = kernels.GaussKernel(tgamma);
            this.SystemKernel = kernels.GaussKernel(sgamma);
            %this.ParamKernel = kernels.GaussKernel(pgamma);
        end
        
        function atd = selectTrainingData(this, modeldata)
            % Selects a subset of the projection training data linearly
            % spaced. The number of samples taken is determined by the
            % ApproxExpansionSize number.
            %
            % Important:
            % Note that the selected training data is projected into the
            % precomputed subspace if spacereduction is performed.
            %
            % Overrides the default method in BaseApprox.
            %
            % See also:
            % models.BaseFullModel.off4_genApproximationTrainData
            
            % Validity checks
            sn = modeldata.ProjTrainData;
            if isempty(sn)
                error('No projection training data available to take approximation training data from.');
            end
            
            selection = round(linspace(1,size(sn,2),...
                    min(this.ApproxExpansionSize,size(sn,2))));
            atd = sn(:,selection);
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
            copy.ApproxExpansionSize = this.ApproxExpansionSize;
        end
    end
    
    methods(Access=protected)
        
        function gen_approximation_data(this, model, xi, ti, mui, fxi)
            % Computes the approximation according to the concrete
            % approximation strategy.
            % Fills the Ma, off and snData properties of the
            % CompwiseKernelCorefun with data.
            
            this.snData.xi = xi;
            this.snData.ti = ti;
            this.snData.mui = mui;
            n = size(xi,2);
            
            %             this.guessKernelConfig;
            %             factors = [.1 .25 .5 .75 1 1.5 2 5 10];
            %             minDiff = Inf;
            %             for idx=1:length(factors)
            %                 gamma = factors(idx)*this.sg;
            %                 this.SystemKernel.Gamma = gamma;
            
            % Call subclass preparation method and pass kernel matrix
            this.prepareApproximationGeneration(...
                this.evaluateAtCenters(xi, ti, mui));
            
            if this.ComputeParallel
                this.computeParallel(model.Data.ApproxfValues);
            else
                this.computeSerial(model.Data.ApproxfValues);
            end
            
            %                 afxi = this.evaluate(xi,ti,mui);
            %                 di = norm(fxi-afxi);
            %                 fprintf('Approx difference on centers for gamma=%f: di\n',gamma);
            %                 if (di < minDiff)
            %                     bestGamma = gamma;
            %                     bestMa = this.Ma;
            %                     bestoff = this.off;
            %                 end
            %             end
            %             this.SystemKernel.Gamma = bestGamma;
            %             this.Ma = bestMa;
            %             this.off = bestoff;
            
            % Reduce the snapshot array and coeff data to the
            % really used ones! This means if any snapshot x_n is
            % not used in any dimension, it is kicked out at this
            % stage.
            hlp = sum(this.Ma ~= 0,1);
            usedidx = find(hlp > 0);
            if length(usedidx) < n
                this.Ma = this.Ma(:,usedidx);
                this.snData.xi = xi(:,usedidx);
                if ~isempty(ti)
                    this.snData.ti = ti(:,usedidx);
                end
                if ~isempty(mui)
                    this.snData.mui = mui(:,usedidx);
                end
            end
            
            % @todo find out when sparse representation is more
            % efficient!
            if sum(hlp) / numel(this.Ma) < .5
                this.Ma = sparse(this.Ma);
            end
            
            % dont use offset vector if none are given
            if all(this.off == 0)
                this.off = [];
            end
            
        end
        
    end
    
    methods(Access=private)
        
        function computeSerial(this, fxi)
            %% Non-parallel execution
            wh = waitbar(0,'Initializing component-wise kernel approximation');
            try
                n = size(this.snData.xi,2);
                fdims = size(fxi,1);
                this.Ma = zeros(fdims, n);
                this.off = zeros(fdims, 1);
                for fdim = 1:fdims
                    waitbar(fdim/fdims+10,wh,sprintf('Computing approximation for dimension %d/%d ... %2.0f %%',fdim,fdims,(fdim/fdims)*100));
                    
                    % Call template method
                    [ai, b, svidx] = this.calcComponentApproximation(fxi(fdim,:));
                    if ~isempty(svidx)
                        this.Ma(fdim,svidx) = ai;
                    else
                        this.Ma(fdim,:) = ai;
                    end
                    this.off(fdim) = b;
                end
                
                close(wh);
            catch ME
                close(wh);
                rethrow(ME);
            end
        end
        
        function computeParallel(this, fxi)
            %% Parallel execution
            n = size(this.snData.xi,2);
            fdims = size(fxi,1);
            fprintf('Starting parallel component-wise approximation computation of %d dimensions on %d workers...\n',fdims,matlabpool('size'));
            parAI = cell(fdims, n);
            parOff = zeros(fdims, 1);
            parSV = cell(1,fdims);
            % Create handle to speedup communications (otherwise the
            % whole object will be copied to each worker)
            cfun = @this.calcComponentApproximation;
            parfor fdim = 1:fdims
                %waitbar(fdim/fdims+10,wh,sprintf('Computing
                %approximation for dimension %d/%d ... %2.0f %%',fdim,fdims,(fdim/fdims)*100));
                % Call template method
                [parAI{fdim}, b, parSV{fdim}] = cfun(fxi(fdim,:));
                parOff(fdim) = b;
            end
            
            this.Ma = zeros(fdims, n);
            for idx = 1:fdims
                if ~isempty(parSV{idx})
                    this.Ma(idx,parSV{idx}) = parAI{idx};
                else
                    this.Ma(idx,:) = parAI{idx};
                end
            end
            this.off = parOff;
        end
    end
    
    methods(Abstract, Access=protected)
        % Preparation template method.
        %
        % Is called once before the component-wise approximation
        % computation. Here subclasses may create the conrete
        % single-dimension approximation algorithm implementations
        %
        % Parameters:
        % K: The Kernel matrix created from the snData `x_i`
        %
        % See also: calcComponentApproximation
        prepareApproximationGeneration(this, K);
        
        % Single dimension approximation computation.
        %
        % Here the concrete class performs the approximation calculation
        % for given function evaluation points `fx_i` at the snData
        % `x_i` also used to compute the kernel matrix for
        % prepareApproximationGeneration. Within this method the class
        % properties cData and csMap are to be
        % manipulated/filled with approximation data. cData ist
        % preallocated to dims x size(xi,2).
        %
        % Parameters:
        % fxi: The function values `f(x_i)` as row vector.
        %
        % Return values:
        % ai: The coefficients `\alpha_{k,i}` of `f_k(x)`.
        % b: The offset for `f_k(x)`.
        % svidx: The used support vector indices `i` of `x_i`. Optional,
        % leave empty if all are used.
        %
        % See also: prepareApproximationGeneration
        [ai, b, svidx] = calcComponentApproximation(this, fxi);
    end
    
end


