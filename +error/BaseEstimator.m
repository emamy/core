classdef BaseEstimator < KerMorObject & ICloneable
    % Base class for all error estimators.
    %
    %   In KerMor any error estimators have the chance to add additionaly
    %   terms to the ODE function by using the
    %   @ref evalODEPart and @ref getE0 functions.
    %   Finally, the method @ref process gets called when the ODE solver
    %   finished to allow for any processing / final computation of error
    %   estimates.
    %
    % @author Daniel Wirtz @date 24.11.2010
    %
    % @change{0,3,sa,2011-04-23} Implemented Setters for the properties of
    % this class other than ReducedModel
    
    properties(Dependent)
        % Flag that indicates whether error estimation is used or not.
        %
        % Default:
        % true
        Enabled;
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
    
    properties(Access=private)
        fEnabled = true;
    end
    
    methods
        
        function setReducedModel(this, rmodel)
            % Performs a validity check for the given model and sets up the
            % estimator for use with the specified model.
            %
            % Override in subclasses and call this method first to perform
            % any additional offline computations
            
            % Calls the static method of the instance (i wonder how long
            % this will work :-)!)
            msg = this.validModelForEstimator(rmodel);
            if ~isempty(msg)
                error(msg);
            end
            % Assign reduced model to local protected property
            this.ReducedModel = rmodel;
        end
        
        function clear(this)
            % Clears the last error set by the estimator.
            this.LastError = [];
        end
        
        function copy = clone(this, copy)
            % Creates a copy of this error estimator.
            % Call only allowed from subclasses.
            if nargin < 2
                error('Clone to this method can only be called from subclasses.');
            end
            copy.ExtraODEDims = this.ExtraODEDims;
            copy.LastError = this.LastError;
            copy.Enabled = this.Enabled;
            % No cloning of the associated reduced model.
            copy.ReducedModel = this.ReducedModel;
        end
              
        function set.LastError(this, value)
            this.LastError = value;
        end
        
        function set.ExtraODEDims(this, value)
            if ~isposintscalar(value) && value ~= 0
                error('The value must be 0 or positive integer');
            end
            this.ExtraODEDims = value;            
        end
        
        function set.ReducedModel(this, value)
            if ~isa(value,'models.ReducedModel')
                error('The given value has to be a models.ReducedModel instance.');
            end
            this.ReducedModel = value;
        end
        
        function set.Enabled(this, value)
            if ~islogical(value)
                error('Enabled property must be a boolean flag');
            elseif ~value
                this.LastError = [];
            end
            this.fEnabled = value;
        end
        
        function value = get.Enabled(this)
            value = this.fEnabled;
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
        
        % Gets the initial condition vector for additional used ODE
        % dimensions.
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
            elseif isempty(error.TPWLLocalLipEstimator.validModelForEstimator(model))
                est = error.TPWLLocalLipEstimator(model);
            elseif isempty(error.GlobalLipKernelEstimator.validModelForEstimator(model))
                est = error.GlobalLipKernelEstimator(model);
            elseif isa(model.System.f,'models.synth.KernelTest')
                est = error.ExperimentalEstimator(model);
            else
                est = error.DefaultEstimator(model);
                fprintf('BaseEstimator::getEstimator: No suitable error estimator found for given model. Using the default estimator (disabled).\n');
            end
        end
        
        function res = test_ErrorEstimators
            % Quick test for estimators.
            demo = RandomModelEstimatorDemo;
            demo.Dims = 3;
            demo.Run;
            res = true;
        end
        
    end
    
    methods(Static, Abstract)
        % Abstract static method that forces subclasses to specify whether
        % an estimator can be used for a given model or not.
        errmsg = validModelForEstimator(rmodel);
    end
end

