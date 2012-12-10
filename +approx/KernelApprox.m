classdef KernelApprox < approx.BaseApprox & dscomponents.ParamTimeKernelCoreFun
% KernelApprox: Base class for component-wise kernel approximations
%
% For each dimension `k` there is a representation
% ``f_k(x) = \sum\limits_{i=1}^N \alpha_{k,i}\Phi(x,x_i)``
% for the approximation. The property Ma contains in row `k` all
% indices `\alpha_k` used in the `k`-th dimension.
%
% The property Centers (inherited from kernels.ParamTimeKernelExpansion) contains all state
% variable Centers that are relevant for the evaluation of the approximated function. (No
% matter how many originally have been used for the approximation computation)
%
% @change{0,5,dw,2011-11-02}
% - Adopted to the new interface for approximation computation (passing data.ApproxTrainData
% instance instead of 'xi,ti,mui' parameters)
% - Also cloning the approximation algorithm instance at KernelApprox.clone
%
% @change{0,5,dw,2011-10-16} Included setting the projection-induced
% norm `V^tGV` to any kernels.ARBFKernel implementing classes
% as required by the theory. (So far made no difference as we had
% `V^tGV=I_r` all the time)
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
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

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
            
            this.Algorithm = approx.algorithms.VKOGA;
             
            this.registerProps('Algorithm');
        end
        
        function approximateSystemFunction(this, model)
            % If V=W, we have W^tV = I_r by assumption, so if G=1 we have
            % V^tGV = I_r and we dont need to set a custom norm for the
            % kernels (would mean additional rounding error)
%             if isa(this.Kernel,'kernels.ARBFKernel') && ...
%                     ~isequal(model.Data.V,model.Data.W) || model.G ~= 1
%                 % Set norm matrix to V^tGV as required by theory.
%                 this.Kernel.G = model.Data.V'*(model.G*model.Data.V);
%             end
            this.Kernel.G = model.G;
                        
            % First argument: this kernel expansion!
            this.Algorithm.computeApproximation(this, model.Data.ApproxTrainData);
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
            
            copy.Algorithm = this.Algorithm.clone;
        end
    end
    
end


