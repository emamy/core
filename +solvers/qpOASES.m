classdef qpOASES < solvers.BaseQPSolver
    % Solves a quadratic program following the active set qp solving
    % algorithms described in the literature:
    %
    % M.J. Best. Applied Mathematics and Parallel Computing, chapter An
    % Algorithm for the Solution of the Parametric Quadratic Programming
    % Problem, pages 57?76. Physica-Verlag, Heidelberg, 1996.
    %
    % H.J. Ferreau. An Online Active Set Strategy for Fast Solution of
    % Parametric Quadratic Programs with Applications to Predictive Engine
    % Control. Master?s thesis, University of Heidelberg, 2006.
    %
    % H.J. Ferreau, H.G. Bock, and M. Diehl. An online active set strategy
    % to overcome the limitations of explicit MPC. International Journal of
    % Robust and Nonlinear Control, 18(8):816?830, 2008.
    %
    % See also: BaseQPSolver
    %
    % @todo write test that ensures correct order of lagrange multipliers
    % over different versions (return variable d; KerMor currently sticks
    % with this order)
    %
    % @author Daniel Wirtz @date 27.10.2010
       
    methods
        
        function this = qpOASES
            % Creates a qpOASES Wrapper for the Matlab-Interface.
            if isempty(which('qpOASES'))
                error('qpOASES installation not found. Added to path?');
            end
            this.Name = 'OASES QP Solver';
        end
        
        function copy = clone(this)
            copy = solvers.qpOASES;
            copy = clone@solvers.BaseQPSolver(this, copy);
        end
        
    end
    methods(Access=protected)
        
        function [p,d,cflag,info] = internal_solve(this,Q,c,lb,ub,A,lbA,ubA,x0)
            
            % Call solver
            if nargin < 9
                % no starting vector given
                [obj,p,d,status,nWSRout] = qpOASES(Q,c,A,lb,ub,lbA,ubA,this.MaxIterations);
            else
                % use provided starting point
                [obj,p,d,status,nWSRout] = qpOASES(Q,c,A,lb,ub,lbA,ubA,this.MaxIterations,x0);
            end
            
            % Transform to KerMor-compliant output
            cflag = status == 0;
            % Append specific outputs to info struct
            info.status = status;
            info.Iterations = nWSRout;
        end
    end
    
end

