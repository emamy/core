classdef DefaultEstimator < error.BaseEstimator
    %DEFAULTESTIMATOR Default error "estimator" for reduced models
    % Standard estimator that is independent from any special reduced model
    % since it computes the full error!
    %
    % @attention This error estimator computes the true '''state space''' error norm.
    % If output is used, the correct state space error norm is mapped to the ''estimated'' output error
    % norm `\no{y} \leq \no{C}\no{e}`, so the computed "correct" output error estimation will
    % '''NOT''' equal the real output error norm `\no{y} = \no{C(x-Vz)}`, as only `\no{e(t)}` but
    % not the vector-valued `e(t)` is known. However, the extra property
    % error.DefaultEstimator.RealOutputError contains the true output error `\no{C(x-Vz)}` as during
    % processing the correct error `e(t)` is known.
    %
    % @author Daniel Wirtz @date 2010-11-24
    %
    % @change{0,4,dw,2011-06-01} Introduced a new property error.DefaultEstimator.RealOutputError to
    % reflect that the actual output error also cannot be obtained by this error estimator. See the
    % class comments for details.
    %
    % @change{0,4,dw,2011-05-23} Adopted to the new error.BaseEstimator interface with separate output
    % error computation.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing    
    
    properties(SetAccess=private)
        % The true output error `\no{C(x(t)-Vz(t))}` for `t\in[0,T]`.
        %
        % In difference to the OutputError property, this field contains the true output error of
        % the reduced simulation. See the class comments for details.
        %
        % See also: OutputError
        RealOutputError;
    end
    
    methods
        function this = DefaultEstimator(rmodel)
            % Creates the default error estimator that works with every
            % model since it computes the full system error.
            %
            % Disabled per default since estimations are expensive.
            this.Enabled = false;
            this.ExtraODEDims = 0;
            
            if nargin == 1
                this.setReducedModel(rmodel);
            end
        end
               
        function copy = clone(this)
            % Clones this DefaultEstimator
            copy = error.DefaultEstimator(this.ReducedModel);
            copy = clone@error.BaseEstimator(this, copy);
            % No local properties there.
        end
        
        function eint = evalODEPart(this, x, t, mu, ut)%#ok
            eint = [];
        end
                
        function e0 = getE0(this, mu)%#ok
            % This error estimator does not use any ODE dimensions, so the
            % initial error part is an empty matrix.
            e0 = [];
        end
        
        function prepareConstants(this, mu, inputidx)%#ok
            % Nothing to do here
        end
    end
    
    methods(Access=protected)
        function postprocess(this, t, x, mu, inputidx)
            m = this.ReducedModel.FullModel;
            % Compute full solution
            [tf,yf,time,xf] = m.simulate(mu,inputidx);
            xr = x(1:end-this.ExtraODEDims,:);
            if ~isempty(this.ReducedModel.V)
                diff = xf-this.ReducedModel.V*xr;
            else
                diff = xf-xr;
            end
            
            % Re-scale
            this.StateError = sqrt(sum(diff.*(this.ReducedModel.GScaled*diff),1));
            
            yr = this.ReducedModel.System.C.computeOutput(t, xr, mu);
            this.RealOutputError = sqrt(sum((yf-yr).^2,1));
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)%#ok
            errmsg = [];
        end
    end
end

