classdef SolverEventData < event.EventData
% SolverEventData: 
%
% Contains information about the actual times and state vectors a solver is about to use or has
% used, depending on the context the event is fired in.
%
% See also: BaseSolver.PreSolve BaseSolver.PostSolve
%
% @author Daniel Wirtz @date 2011-05-31
%
% @new{0,4,dw,2011-05-31} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        Times;
        States = [];
    end
    
    methods
        function this = SolverEventData(t,x)
            if nargin > 0
                this.Times = t;
                if nargin > 1
                    this.States = x;
                end
            end
        end
    end
    
end