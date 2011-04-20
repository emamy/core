classdef IJacobian < handle
% IJacobian: Interface for all core functions that provide an analytical jacobian matrix for the
% state components
%
% Used for supplying a jacobian matrix in implicit solvers, for example. One nice thing about kernel
% expansions is the fact that an jacobian matrix is readily computable, so any
% dscomponents.CompwiseKernelCoreFun natively implements this interface.
%
% @author Daniel Wirtz @date 2011-04-15
%
% @new{0,3,dw,2011-04-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
%     properties(Dependent)
%         hStateJacobian;
%     end
%     
%     methods
%         function h = get.hStateJacobian(this)
%             h = @this.getStateJacobian;
%         end
%     end
        
    methods(Abstract)
        % Returns the jacobian matrix with respect to the state variable at the given point.
        getStateJacobian(x,t,mu);
    end
    
end