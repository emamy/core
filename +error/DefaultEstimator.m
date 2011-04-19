classdef DefaultEstimator < error.BaseEstimator
    %DEFAULTESTIMATOR Default error "estimator" for reduced models
    % Standard estimator that is independent from any special reduced model
    % since it computes the full error!
    %
    % Warning:
    % If a nontrivial output conversion is used (i.e. `C\neq1`), the
    % default estimator computes `e_y(t) = ||C||\cdot||e(t)||`, which will
    % result in an higher error estimate than the true error.
    
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
        
        function process(this, t, x, mu, inputidx)
            m = this.ReducedModel.FullModel;
            % Compute full solution
            [tf,xf] = m.computeTrajectory(mu,inputidx);
            xr = x(1:end-this.ExtraODEDims,:);
            if ~isempty(this.ReducedModel.V)
                diff = xf-this.ReducedModel.V*xr;
            else
                diff = xf-xr;
            end
            
            % Re-scale
            %diff = bsxfun(@times,diff,this.ReducedModel.System.StateScaling);
            y = sqrt(sum(diff.*(this.ReducedModel.GScaled*diff),1));
            
            % Convert to exact error on output level
%             diffy = m.System.C.computeOutput(t,diff,mu);
%             this.LastError = sqrt(sum(diffy.*(this.ReducedModel.G*diffy),1));
%             return;
            
            % Get error
            if m.System.C.TimeDependent
                for idx=1:length(t)
                    y(idx) = norm(m.System.C.evaluate(t(idx),mu))*y(idx);
                end
            else
                y = norm(m.System.C.evaluate([],mu))*y;
            end
            this.LastError = y;
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

