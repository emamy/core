classdef Base < handle
% Base: Interface for initial error computing classes.
%
% @author Daniel Wirtz @date 2011-07-04
%
% @new{0,5,dw,2011-07-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods(Abstract)
        % Returns the initial value error for given parameter `\mu`
        %
        % Template method.
        e0 = getE0(this, mu);
        
        % Sets the effectively used reduced model.
        prepareForReducedModel(this, rm);
    end
    
end