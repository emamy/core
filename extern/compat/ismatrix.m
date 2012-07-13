function ismat = ismatrix(value)
% ismatrix: Compatibility function for matlab versions smaller than 2011b
%
%
%
% @author Daniel Wirtz @date 2012-04-06
%
% @new{0,6,dw,2012-04-06} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
if exist('ismatrix','builtin') == 5
    ismat = builtin('ismatrix', value);
else
    ismat = isscalar(value) || ndims(value) == 2 && length(value) > 1;
end
end