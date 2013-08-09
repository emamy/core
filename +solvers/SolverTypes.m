classdef SolverTypes
% SolverTypes: Enumeration class classifying a ode solver type as explicit solver, implicit solver or Matlab solver.
%
% A maximal integration step is only relevant for exlicit solvers.
%
% @new{0,7,ts,2013-07-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    enumeration
        Explicit, Implicit, MLSolver
    end
end