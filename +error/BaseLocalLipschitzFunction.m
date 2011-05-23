classdef BaseLocalLipschitzFunction < KerMorObject
% BaseLocalLipschitzFunction: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=protected)
        bellfcn;
    end
    
    methods
        
        function this = BaseLocalLipschitzFunction(bellfunc)
            this = this@KerMorObject;
            this.bellfcn = bellfunc;
        end
        
    end
    
    methods(Abstract)
        ci = evaluate(this,di,Ct,t,mu);
    end
    
end