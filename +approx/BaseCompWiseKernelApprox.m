classdef BaseCompWiseKernelApprox < approx.BaseApprox & ...
        dscomponents.CompwiseKernelCoreFun & IParallelizable & ...
        approx.IAutoConfig
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
    % See also: BaseKernelApprox IAutoConfig

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
        
        function autoconfigure(this, model)
            % Implements the template method from IAutoConfigure.
            %
            % For this class autoconfiguration means detection of the
            % "ideal" radius for gaussian kernels, if used. The strategy is
            % to enforce that for the largest distance between any two
            % considered centers the sum of both nearby kernel evaluations
            % equals one, i.e.
            % ``e^{-\frac{\left(\frac{d}{2}\right)^2}{\gamma}} =
            % \frac{1}{2}``
            % if `d` is the largest distance.
            % 
            % Parameters:
            % model: The current model instance
            %
            % See also: IAutoConfigure
            
            % Settings.
            %zero = 1e-4;
            %trange = 3; % nonzero over trange times the dt-distance
            %srange = 3; % nonzero over srange times the maximum distance within the training data
            %prange = 2; % nonzero over prange times the param samples distance
            data = model.Data.ApproxTrainData;
            v = unique(round(data(1,:)));
            
            %% State kernel gamma
            if isa(this.SystemKernel,'kernels.GaussKernel')
                xd = sqrt(sum(data(4:end,:).^2));

                % Find samples for each parameter    
                maxdiff = zeros(1,length(v));
                for muidx = 1:length(v)
                    sel = data(1,:) == v(muidx);
                    tmp = xd(sel);
                    maxdiff(muidx) = max(abs(tmp(1:end-1)-tmp(2:end)));
                end
                d = max(maxdiff);
                this.SystemKernel.setGammaForDistance(d/2,.5);
            end
            
            %% Time kernel
            if isa(this.TimeKernel,'kernels.GaussKernel')
                warning('Code:unchecked','Implementation not yet finished/ideal!');
                this.TimeKernel.setGammaForDistance(model.dt/2,.5);
            end
            
            %% Param kernel
            if isa(this.ParamKernel,'kernels.GaussKernel')
                params = model.Data.getParams(v);
                % Gives a matrix with parameters in each column
                mud = sqrt(sum(params.^2));
                % Create distance matrix
                dist = abs(repmat(mud,size(mud,2),1)-repmat(mud',1,size(mud,2)));
                this.ParamKernel.setGammaForDistance(max(dist(:))/2,.5);
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


