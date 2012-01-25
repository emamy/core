classdef IDemoProvider < handle
% IDemoProvider: Interface for any KerMor class that provides an executable standalone-demo.
%
% Implementing classes just have to call runDemo and can go right ahead!
%
% @author Daniel Wirtz @date 2012-01-25
%
% @new{0,6,dw,2012-01-25} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Abstract, Static)
        % Demo execution main function.
        runDemo;
    end
    
end