classdef StopFlag
% StopFlag: Flags that algorithms can use to specify the reason for their termination.
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-02-20
%
% @new{0,7,dw,2013-02-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Constant)
        % Absolute error tolerance reached
        ABS_ERROR = 1;
        
        % Relative error tolerance reached
        REL_ERROR = 2;
        
        % Maximum size reached
        MAX_SIZE = 3;
        
        % Stop flag from the VKOGA algorithm.
        NEGATIVE_POWFUN = 4;
        
        % Maximum number of iterations reached
        MAX_ITER = 5;
        
        % Algorithm terminated otherwisely successful
        SUCCESS = 6;
        
        % SVR-specific flag that indicates that no support vectors have been found.
        NO_SUPPORT_VECTORS_FOUND = 7;
        
        % SVR-SMO specific flag. Indicates that the E+T values are smaller than the prescribed
        % tolerance.
        TOL_OK = 8;
    end
  
end