classdef DEIM < approx.BaseApprox & general.DEIM
% DEIM: Wrapper for KerMor dynamical systems of the general.DEIM class
%
% @author Daniel Wirtz @date 2012-03-26
%
% @new{0,6,dw,2012-05-30} Added a custom implementation of the
% getStateJacobian method
%
% @change{0,6,dw,2012-05-29} Besides many work-in-progress changes, now the
% evaluate function is implemented directly. further, some matrices needed
% for the DEIMEstimator have been included and are being computed.
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% @todo think of exporting the jrow, jend, S properties to ACompEvalCoreFun
% and precompute stuff there; projection happens at that stage
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = DEIM
            this = this@general.DEIM;
            this = this@approx.BaseApprox;
            
            this.CustomProjection = true;
            this.TimeDependent = true;
        end
        
        function approximateSystemFunction(this, model)
            this.computeDEIM(model.System.f, model.Data.ApproxTrainData.fxi);
        end
        
        function fx = evaluate(this, x, t, mu)
            fx = evaluate@general.DEIM(this, x, t, mu);
        end
        
        function J = getStateJacobian(this, x, t, mu)
            J = getStateJacobian@general.DEIM(this, x, t, mu);
            % Finite difference one:
            %J2 = getStateJacobian@dscomponents.ACoreFun(this, x, t, mu);
        end
        
        function fx = evaluateCoreFun(varargin)
            % do nothing as evaluate is overridden directly
            error('Do not call me!');
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected = project@approx.BaseApprox(this, V, W, projected);
            projected = project@general.DEIM(this, V, W, projected);
        end
        
        function copy = clone(this)
            copy = approx.DEIM;
            copy = clone@approx.BaseApprox(this, copy);
            copy = clone@general.DEIM(this, copy);
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, varargin)
            obj = loadobj@general.DEIM(obj, varargin{:});
            obj = loadobj@approx.BaseApprox(obj, varargin{:});
        end
    end
end