classdef BaseEstimator < KerMorObject & ICloneable & ISimConstants
    % Base class for all error estimators.
    %
    %   In KerMor any error estimators have the chance to add additionaly
    %   terms to the ODE function by using the
    %   @ref evalODEPart and @ref getE0 functions.
    %   Finally, the method @ref postProcess gets called when the ODE solver
    %   finished to allow for any processing / final computation of error
    %   estimates.
    %
    % @author Daniel Wirtz @date 24.11.2010
    %
    % @change{0,5,dw,2011-07-07} Added a new property e0Comp that computes the initial error for the
    % given reduced model, depending on its initial value class. In consequence, the getE0 method is
    % now found at this base error estimator class.
    %
    % @change{0,4,dw,2011-05-23} Created a new interface with separate output
    % error computation. Postprocessing is now a template method and a new property
    % error.BaseEstimator.OutputError has been introduced. Renamed the LastError property to
    % error.BaseEstimator.StateError.
    %
    % @change{0,3,sa,2011-04-23} Implemented Setters for the properties
    %
    
    properties(Dependent)
        % Flag that indicates whether error estimation is used or not.
        %
        % @default true
        Enabled;
    end
    
    properties(SetAccess=protected)
        % The reduction state-space error from the last simulation
        %
        % @default []
        %
        % See also: models.ReducedModel#ErrorEstimator
        StateError = [];
               
        % The dimensions added to the ODE function by the estimator.
        %
        % Please set to correct value in subclasses as simulations are not
        % possible if set incorrect.
        %
        % @default 0
        ExtraODEDims = 0;
    end
    
    properties(SetAccess=private)
        % The output error from the last simulation
        %
        % @default []
        %
        % See also: models.ReducedModel#ErrorEstimator
        OutputError = [];
    end
    
    properties(Access=protected)
        % The reduced model associated with the error estimator.
        %
        % @type models.ReducedModel
        ReducedModel;
    end
    
    properties(Access=private)
        fEnabled = true;
        e0Comp;
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
                builtin('error',msg);
            end
            % Assign reduced model to local protected property
            this.ReducedModel = rmodel;
            
            if isa(rmodel.System.x0,'dscomponents.AffineInitialValue')
                this.e0Comp = error.initial.AffineParametric(rmodel);
            else
                this.e0Comp = error.initial.Constant(rmodel);
            end
        end
        
        function postProcess(this, x, t, mu, inputidx)
            % Post-processes the error estimator ODE part after reduced simulation computation
            %
            % Here the OutputError fields
            %
            % Parameters:
            % t: The times at which the reduced simulation was computed
            % x: The reduced simulation's system state PLUS the error
            % estimation values in the last this.ExtraODEDims rows.
            % mu: The current parameter `\mu`
            % inputidx: The current input index
            
            % Call template method to compute the state error in subclasses
            this.postprocess(x, t, mu, inputidx);
            
            if all(this.StateError == 0)
                warning('BaseEstimator:postProcess','State error is all zero. Attention!');
            end
            
            % Tranform to output error estimation (if used)
            C = this.ReducedModel.FullModel.System.C;
            if ~isempty(C)
                % Get error
                if C.TimeDependent
                    e = this.StateError;
                    for idx=1:length(t)
                        e(idx) = norm(m.System.C.evaluate(t(idx),mu))*e(idx);
                    end
                else
                    this.OutputError = norm(C.evaluate([],mu))*this.StateError;
                end
            else
                this.OutputError = this.StateError;
            end
        end
        
        function clear(this)
            % Clears the last error set by the estimator.
            this.StateError = [];
            this.OutputError = [];
        end
        
        function copy = clone(this, copy)
            % Creates a copy of this error estimator.
            % Call only allowed from subclasses.
            if nargin < 2
                error('Clone to this method can only be called from subclasses.');
            end
            copy.ExtraODEDims = this.ExtraODEDims;
            copy.StateError = this.StateError;
            copy.OutputError = this.OutputError;
            copy.Enabled = this.Enabled;
            % No cloning of the associated reduced model.
            copy.ReducedModel = this.ReducedModel;
            copy.e0Comp = this.e0Comp;
        end
        
        function e0 = getE0(this, mu)
            % Calls the inner initial error computation strategy.
            e0 = this.e0Comp.getE0(mu);
        end
        
    end
    
    %% Getter & Setter
    methods
              
        function set.StateError(this, value)
            this.StateError = value;
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
                this.StateError = [];
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
    end
    
    methods(Abstract, Access=protected)
        % Allows for any post-processing after the ODE solver finished
        % integration of the system/error system part.
        %
        % Parameters:
        % t: The times at which the reduced simulation was computed
        % x: The reduced simulation's system state PLUS the error
        % estimation values in the last this.ExtraODEDims rows.
        postprocess(this, t, x, mu, inputidx);
    end
    
    methods(Static)
        function est = getEstimator(model)
            % Factory method that creates a suitable error estimator for
            % the given model.
            %
            % Tries to always select the best estimator available for the
            % model. Of course the error estimator can be changed manually
            % later on.
            %
            % @todo overhaul this method of assigning an error estimator to a reduced model.
            % at this method local knowledge of all available error estimators has to be present
            % anyways if going through this in a if then fashion. otherwise, see if reflection may
            % be used here!
            if isempty(error.IterationCompLemmaEstimator.validModelForEstimator(model))
                est = error.IterationCompLemmaEstimator(model);
%             elseif isempty(error.TPWLLocalLipEstimator.validModelForEstimator(model))
%                 est = error.TPWLLocalLipEstimator(model);
            elseif isempty(error.GLEstimator.validModelForEstimator(model))
                est = error.GLEstimator(model);
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
            demo.start;
            res = true;
        end
        
    end
    
    methods(Static, Abstract)
        % Abstract static method that forces subclasses to specify whether
        % an estimator can be used for a given model or not.
        errmsg = validModelForEstimator(rmodel);
    end
end

