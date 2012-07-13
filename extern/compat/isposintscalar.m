function res = isposintscalar(value)
% isposintscalar: Backwards-compatibility function for matlab versions greater than 2012a
%
% @author Daniel Wirtz @date 2012-07-03
%
% @new{0,6,dw,2012-07-03} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    if exist('isposintscalar','builtin') == 5
        res = builtin('isposintscalar', value);
    else
        res = isposrealscalar(value) && (round(value) - value) < sqrt(eps);
    end
end