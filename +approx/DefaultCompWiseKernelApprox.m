classdef DefaultCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
    % Default component-wise kernel approximation algorithm
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
        function computeCompwiseApprox(this, model, fx)
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
            
            
            % Call coeffcomp preparation method and pass kernel matrix
            this.CoeffComp.init(this.getKernelMatrix);
            
            % Call protected method
            this.computeCoeffs(fx);
             
            % dont use offset vector if none are given
            if all(this.off == 0)
                this.off = [];
            end     
        end          
    end
end


