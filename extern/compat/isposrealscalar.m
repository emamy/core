function res = isposrealscalar(value)
% isposintscalar: Backwards-compatibility function for matlab versions greater than 2012a
%
% @author Daniel Wirtz @date 2012-07-03
%
% @new{0,6,dw,2012-07-03} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    if exist('isposrealscalar','builtin') == 5
        res = builtin('isposrealscalar', value);
    else
        res = isscalar(value) && value >= 0;
    end
end