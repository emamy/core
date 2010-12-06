classdef DefaultEstimator < error.BaseEstimator
    %DEFAULTESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
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
        
        function offlineComputations(this)%#ok
            % nothing to do here
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
        
        function process(this, t, x, mu, inputidx)%#ok
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

