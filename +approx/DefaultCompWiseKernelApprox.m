classdef DefaultCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
    % Default component-wise kernel approximation algorithm
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % @change(0,3,sa,2011-04-21) Implemented Setter for the class property
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
        
        function approximateCoreFun(this, model)
            % Load snapshots
            atd = model.Data.ApproxTrainData;
            
            % Compile necessary data
            xi = atd(4:end,:);
            ti = atd(3,:);
            muidx = atd(1,:);
            if all(muidx == 0)
                mui = [];
            else
                mui = model.Data.ParamSamples(:,muidx);
            end
            % Set AKernelCoreFun centers
            this.Centers.xi = xi;
            this.Centers.ti = ti;
            this.Centers.mui = mui;
            n = size(xi,2);
            
            % Call coeffcomp preparation method and pass kernel matrix
            this.CoeffComp.init(this.getKernelMatrix);
            
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
        
        function set.ApproxExpansionSize(this, value)
            if ~isposintscalar(value)
                error('The value should be a positive integer');
            end
            this.ApproxExpansionSize = value;
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
end


