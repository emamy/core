classdef BaseEstimator < handle
    % Base class for all error estimators.
    %
    %   In KerMor any error estimators have the chance to add additionaly
    %   terms to the ODE function by using the
    %   @ref evalODEPart and @ref getE0 functions.
    %   Finally, the method @ref process gets called when the ODE solver
    %   finished to allow for any processing / final computation of error
    %   estimates.
    
    properties
        % Flag that indicates whether error estimation is used or not.
        %
        % Default:
        % true
        Enabled = true;
    end
    
    properties(SetAccess=protected)
        % The reduction error data from the last simulation
        %
        % Default:
        % []
        %
        % See also: models.ReducedModel#ErrorComputation
        LastError = [];
        
        % The dimensions added to the ODE function by the estimator.
        %
        % Please set to correct value in subclasses as simulations are not
        % possible if set incorrect.
        %
        % Default:
        % 0
        ExtraODEDims = 0;
    end
    
    properties(Access=protected)
        % The reduced model associated with the error estimator.
        ReducedModel;
    end
    
    methods
        function this = BaseEstimator(rmodel)
            this.ReducedModel = rmodel;
        end
        
        function clear(this)
            % Clears the last error set by the estimator.
            this.LastError = [];
        end
    end
    
    methods(Abstract)
        % Parameters:
        % x: The full extended state variable vector. Extended means that
        % the last @ref ExtraODEDims rows contain the error estimators own
        % data. If not used, implementers must take care to ditch those
        % values if any function evaluations are performed within the
        % integral part.
        % ut: The value of the input function `u(t)` if given, [] else.
        eint = evalODEPart(this, x, t, mu, ut);
        
        % Allows for any post-processing after the ODE solver finished
        % integration of the system/error system part.
        %
        % Parameters:
        % t: The times at which the reduced simulation was computed
        % x: The reduced simulation's system state PLUS the error
        % estimation values in the last this.ExtraODEDims rows.
        process(this, t, x, mu, inputidx);
        
        % Gets the initial vector for additionally used ODE dimensions.
        e0 = getE0(this, mu);
    end
    
    methods(Static)
        function est = getEstimator(model)
            % Factory method that creates a suitable error estimator for
            % the given model.
            %
            % Tries to always select the best estimator available for the
            % model. Of course the error estimator can be changed manually
            % later on.
            if isempty(error.LocalLipKernelEstimator.validModelForEstimator(model))
                est = error.LocalLipKernelEstimator(model);
            elseif isempty(error.GlobalLipKernelEstimator.validModelForEstimator(model))
                est = error.GlobalLipKernelEstimator(model);
            elseif isa(model.System.f,'models.synth.KernelTest')
                est = error.ExperimentalEstimator(model);
            else
                est = error.DefaultEstimator(model);
                warning('error:BaseEstimator','No suitable error estimator found for given model. Using the default estimator (disabled).');
            end
        end
        
        function res = test_ErrorEstimators
            % Quick test for estimators.
            %demo = LocalLipEstimatorDemo(10);
            %demo.Run;
            res = true;
        end
        
    end
    
    methods(Static, Abstract)
        % Abstract static method that forces subclasses to specify whether
        % an estimator can be used for a given model or not.
        errmsg = validModelForEstimator(rmodel);
    end
    
    %% Save & Load
%     methods
%         function s = saveobj(obj)
%             s.Enabled = obj.Enabled;
%             
%         end
%     end
    
%     methods (Static)
%         function obj = loadobj(s)
%             obj.Enabled = s.Enabled;
%         end
%     end
    
end

