classdef SolverTypes
% SolverTypes: Enumeration class classifying a ode solver type as explicit solver, implicit solver or Matlab solver.
%
% A maximal integration step is only relevant for exlicit solvers.
%
% @new{0,7,ts,2013-07-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    enumeration
        Explicit,
        Implicit,
        Unknown
    end
end