classdef KernelApprox < approx.BaseApprox & ...
        dscomponents.ParamTimeKernelCoreFun %& approx.IAutoConfig
    %Base class for component-wise kernel approximations
    %
    % For each dimension `k` there is a representation
    % ``f_k(x) = \sum\limits_{i=1}^N \alpha_{k,i}\Phi(x,x_i)``
    % for the approximation. The property Ma contains in row `k` all
    % indices `\alpha_k` used in the `k`-th dimension.
    %
    % The property Centers (inherited from
    % kernels.ParamTimeKernelExpansion) contains all state variable
    % Centers that are relevant for the evaluation of the approximated
    % function. (No matter how many originally have been used for the
    % approximation computation)
    %
    % @change{0,5,dw,2011-07-07} Now inherits from kernels.ParamTimeKernelExpansion and has a
    % strategy pattern class reference for the approximation generation algorithm "Algorithm".
    % This has been changed so that standalone use of approximation algorithms can also be applied;
    % now this class simply wraps this functionality in the context of the original model reduction
    % scheme.
    %
    % @change{0,4,dw,2011-05-19} Disconnected the Approx classes from taking a BaseModel instance at
    % approx computation. This way external tools can use the approximation algorithms, too.
    %
    % @change{0,4,dw,2011-05-03} Removed off-property dependencies in computeCoeffs routines as
    % expansion offsets are no longer used.
    %
    % @change{0,3,dw,2011-04-26} Implementing the approximateData in this class so that rescaling
    % of the function f values `f_{x_i}` can be used. Subclasses now implement the template method
    % approx.KernelApprox.computeCompwiseApprox.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    %
    % @change{0,3,dw,2011-03-31} 
    % - Added a new property approx.KernelApprox.CoeffComp
    % which allows to use the strategy pattern for computing kernel
    % expansion coefficients.
    % - Moved some old methods into the new class
    % approx.algorithms.DefaultCompWiseKernelApprox to retain the old functionality
    % which selects a subset of the projection training data by linspace
    % selection as approximation data.
    % - Changed method names for computeSerial/Parallel to
    % computeCoeffsSerial/computeCoeffsParallel and visibility to protected
    % as they can be used by any subclass.
    %
    % See also: KernelApprox

    properties(SetObservable)        
        % The algorithm used to create the kernel expansion.
        %
        % @type approx.algorithms.Algorithm
        %
        % @propclass{critical} The approximation strategy for kernel expansions is essential.
        Algorithm;
    end
    
    methods
        
        function this = KernelApprox
            this = this@approx.BaseApprox;
            this = this@dscomponents.ParamTimeKernelCoreFun;
            
            this.Algorithm = approx.algorithms.AdaptiveCompWiseKernelApprox;
             
            this.registerProps('Algorithm');
        end
        
        function approximateSystemFunction(this, model)
            atd = model.Data.ApproxTrainData;
                        
            % First argument: this kernel expansion!
            this.Algorithm.computeApproximation(this, atd.xi, atd.ti, atd.mui, atd.fxi);
        end
                
        function copy = clone(this)
            % Clones the instance.
            %
            % Note:
            % Since we use multiple inheritance we have to call the clone
            % method from both above superclasses. In this case this leads
            % to a double execution of the clone method from
            % dscomponents.ACoreFcn, but this is rather done than ommitted
            % and causing trouble later on.
            copy = approx.KernelApprox;
                       
            copy = clone@dscomponents.ParamTimeKernelCoreFun(this, copy);
            copy = clone@approx.BaseApprox(this, copy);
            
            copy.Algorithm = this.Algorithm;
        end
    end
    
end


