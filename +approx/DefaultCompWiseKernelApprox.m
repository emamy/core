classdef DefaultCompWiseKernelApprox < approx.BaseCompWiseKernelApprox
    % Default component-wise kernel approximation algorithm
    %
    % @author Daniel Wirtz @date 2011-03-31
    %
    % @change{0,4,dw,2011-05-20} Removed the ApproxExpansionSize property as this is now determined
    % by the approx.TrainDataSelector property
    %
    % @change(0,3,sa,2011-04-21) Implemented Setter for the class property
    %
    % See also: BaseApprox BaseCompWiseKernelApprox
    %
    % @new{0,3,dw,2011-03-31} Added this class to keep old approximation
    % generation method.
    
%     properties
%         % The number of projection training data snapshots used to compile
%         % the approximation training data set. So far, the default strategy
%         % implemented in this class simply uses linspace to select a subset
%         % of the specified size.
%         %
%         % Default: 120
%         ApproxExpansionSize = 120;
%     end
    
    methods
%         function set.ApproxExpansionSize(this, value)
%             if ~isposintscalar(value)
%                 error('The value should be a positive integer');
%             end
%             this.ApproxExpansionSize = value;
%         end
        
        function target = clone(this)
            % Clones the instance.
            
            % Create instance as this is the final class so far. If
            % subclassed, this clone method has to be given an additional
            % target argument.
            target = approx.DefaultCompWiseKernelApprox;
            
            target = clone@approx.BaseCompWiseKernelApprox(this, target);
%             % copy local props
%             copy.ApproxExpansionSize = this.ApproxExpansionSize;
        end
    end

    methods(Access=protected) 
        function computeCompwiseApprox(this, xi, ti, mui, fxi)
            % Set AKernelCoreFun centers
            this.Centers.xi = xi;
            this.Centers.ti = ti;
            this.Centers.mui = mui;
            
            % Call coeffcomp preparation method and pass kernel matrix
            this.CoeffComp.init(this.getKernelMatrix);
            
            % Call protected method
            this.computeCoeffs(fxi);  
        end          
    end
end


