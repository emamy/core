classdef DefaultCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
    % Default component-wise kernel approximation algorithm
    %
    % For each dimension `k` there is a representation
    % ``f_k(x) = \sum\limits_{i=1}^N \alpha_{k,i}\Phi(x,x_i) + b_k``
    % for the approximation. The property cData contains in row `k` all
    % indices `\alpha_k` used in the `k`-th dimension. off contains all
    % `b_k`.
    %
    % The property Centers contains all state variable Centers that
    % are relevant for the evaluation of the approximated function. (No
    % matter how many originally have been used for the approximation
    % computation!)
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % See also: BaseApprox BaseCompWiseKernelApprox
    %
    % @new{0,3,dw,2011-03-31} Added this class to keep old approximation
    % generation method.
    
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
        
        function target = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            target = approx.DefaultCompWiseKernelApprox;
            
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
            % copy local props
            copy.ApproxExpansionSize = this.ApproxExpansionSize;
        end
       
    end
    
    methods(Access=protected)
        
        function gen_approximation_data(this, model, xi, ti, mui)
            % Computes the approximation according to the concrete
            % approximation strategy.
            % Fills the Ma, off and Centers properties of the
            % CompwiseKernelCorefun with data.
            
            this.Centers.xi = xi;
            this.Centers.ti = ti;
            this.Centers.mui = mui;
            n = size(xi,2);
            
            %             this.guessKernelConfig;
            %             factors = [.1 .25 .5 .75 1 1.5 2 5 10];
            %             minDiff = Inf;
            %             for idx=1:length(factors)
            %                 gamma = factors(idx)*this.sg;
            %                 this.SystemKernel.Gamma = gamma;
            
            % Call coeffcomp preparation method and pass kernel matrix
            this.CoeffComp.init(...
                this.evaluateAtCenters(xi, ti, mui));
            
            % Call protected method
            this.computeCoeffs(model.Data.ApproxfValues);
            
            % Reduce the snapshot array and coeff data to the
            % really used ones! This means if any snapshot x_n is
            % not used in any dimension, it is kicked out at this
            % stage.
            hlp = sum(this.Ma ~= 0,1);
            usedidx = find(hlp > 0);
            if length(usedidx) < n
                this.Ma = this.Ma(:,usedidx);
                this.Centers.xi = xi(:,usedidx);
                if ~isempty(ti)
                    this.Centers.ti = ti(:,usedidx);
                end
                if ~isempty(mui)
                    this.Centers.mui = mui(:,usedidx);
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
end


