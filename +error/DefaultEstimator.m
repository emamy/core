classdef DefaultEstimator < error.BaseEstimator
    %DEFAULTESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function this = DefaultEstimator(rmodel)
            % Disable per default (expensive!)
            this = this@error.BaseEstimator(rmodel);
            this.Enabled = false;
            this.ExtraODEDims = 0;
        end
        
        function eint = evalODEPart(this, x, t, mu, ut)%#ok
            eint = [];
        end
        
        function process(this, t, x, mu, inputidx)
            m = this.ReducedModel.FullModel;
            [tf,xf] = m.computeTrajectory(mu,inputidx);
            xr = x(1:end-this.ExtraODEDims,:);
            diff = xf- this.ReducedModel.V*xr;
            this.LastError = sqrt(sum(diff.^2));
        end
        
        function e0 = getE0(this, mu)%#ok
            % This error estimator does not use any ODE dimensions, so the
            % initial error part is an empty matrix.
            e0 = [];
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)%#ok
            errmsg = [];
        end
    end
end

