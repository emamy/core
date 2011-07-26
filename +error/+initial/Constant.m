classdef Constant < error.initial.Base
% Constant: 
%
%
%
% @author Daniel Wirtz @date 2011-07-04
%
% @new{0,5,dw,2011-07-04} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        e0;
    end
    
    methods
        function this = Constant(rmodel)
            x0 = rmodel.FullModel.System.x0.evaluate([]);
            % Only project if projection is used
            if ~isempty(rmodel.V)
                x0 = x0 - rmodel.V*(rmodel.W'*x0);
            end
            this.e0 = sqrt(x0'*rmodel.GScaled*x0);
        end
        
        function e0 = getE0(this, mu)%#ok
            e0 = this.e0;
        end
    end
    
end