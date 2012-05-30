classdef AImplSolver < KerMorObject
% AImplSolver: Base abstract class for fully implicit ode solvers
%
% Some models like the models.pcd.PCDModel have a small CFL constant due to fast diffusion. Thus,
% implicit solvers are required to avoid long computation times.
%
% To not have to change the solve interface already present for explicit solvers, this class
% implements the solve method itself and uses a template method to perform the actual implicit
% solving.
%
% @author Daniel Wirtz @date 2011-04-19
%
% @change{0,6,dw,2011-12-07} Renamed this class from BaseImplSolver to
% AImplSolver, now merely providing jacobian matrix information and a
% checkable interface for implicit solvers
%
% @new{0,6,dw,2011-11-27} New property AImplSolver.JPattern that allows to set a sparsity
% pattern for the ode function.
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @new{0,3,dw,2011-04-19} Added this class to have a base for implicit solvers and wrap to the
% normal odefun interface from the explicit ones.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % A function handle to compute the core function's jacobian
        %
        % Optional, any implementations must also work with this property set to [].
        %
        % The function handles parameters must be:
        % - \c t The current time `t`
        % - \c x The current state space vector for time `t`
        % Returned must be the jacobian matrix `\frac{\partial f}{\partial x}(t,x)` at the point
        % `x(t)` and time `t`.
        %
        % @propclass{important} If available, supply a jacobian function evaluation handle to
        % improve speed and reliability of implicit solvers.
        %
        % @type function_handle @default []
        %
        % See also: odeset
        JacFun = [];
        
        % The sparsity pattern of the jacobian `\nabla_x f(x,t,\mu)`
        %
        % @propclass{important} Providing a sparsity pattern to implicit solvers might be
        % crucial for the performance of the solver due to memory restrictions.
        %
        % @type sparsematrix @default []
        JPattern;
    end
    
%     properties(Access=private)
%         fJP = [];
%         fJPD = [];
%     end
    
    methods
        
        function this = AImplSolver
            this.registerProps('JacFun','JPattern');
        end
        
        function set.JacFun(this, value)
            % Sets the jacobian function handle.
            %
            % Parameters:
            % value: a function handle that takes two arguments, `t` and `x`
            if ~isempty(value) && ~isa(value,'function_handle')
                error('JacFun must be a function handle');
            elseif ~isempty(value) && nargin(value) ~= 2
                error('JacFun must take exactly two arguments: t,x');
            end
            this.JacFun = value;
        end
        
        function set.JPattern(this, value)
            if ~isempty(value) && ~issparse(value)
                error('JPattern must be a sparse matrix.');
            end
            this.JPattern = value;
        end
    end
end