classdef BaseEstimator < ICloneable
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
    
    properties
        % Flag that indicates whether error estimation is used or not.
        %
        % Default:
        % true
        Enabled = true;
    end
    
    properties(SetAccess=private, GetAccess=protected)
        M1 = [];
        M2 = [];
        M3 = [];
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
        
        function setReducedModel(this, rmodel)
            % Performs a validity check for the given model and sets up the
            % estimator for use with the specified model.
            
            % Calls the static method of the instance (i wonder how long
            % this will work :-)!)
            msg = this.validModelForEstimator(rmodel);
            if ~isempty(msg)
                error(msg);
            end
            % Assign reduced model to local protected property
            this.ReducedModel = rmodel;
            
            % Perform any offline computations/preparations
            % Only prepare matrices if projection is used
            if ~isempty(this.ReducedModel.V) && ~isempty(this.ReducedModel.W)
                
                % Obtain the correct snapshots
                % Standard case: the approx function is a kernel expansion. it
                % can also be that the system's core function is already a
                % kernel expansion
                fm = this.ReducedModel.FullModel;
                if ~isempty(fm.Approx)
                    % Get full d x N coeff matrix of approx function
                    Ma = fm.Approx.Ma;
                else
                    % Get full d x N coeff matrix of core function
                    Ma = fm.System.f.Ma;
                end
                
                % Compute projection part matrices, without creating a
                % d x d matrix (too big!)
                M = this.ReducedModel.V*(this.ReducedModel.W'*Ma);
                G1 = Ma'*this.ReducedModel.G;
                this.M1 = G1*Ma - 2*G1*M + M'*(this.ReducedModel.G*M);
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(fm.System.B)
                    try
                        B = fm.System.B.evaluate([],[]);
                    catch ME%#ok
                        B = fm.System.B.evaluate(0,this.ReducedModel.System.getRandomParam);
                        warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                    end
                    
                    B2 = this.ReducedModel.V*(this.ReducedModel.W'*B);
                    G2 = B'*this.ReducedModel.G;
                    this.M2 = 2*(G1*B - M'*G2' - G1*B2 + M'*(this.ReducedModel.G*B2));
                    this.M3 = G2*B - 2*G2*B2 + B2'*(this.ReducedModel.G*B2);
                    clear B2 G2;
                end
                clear M G1;
            else
                % No projection means no projection error!
                this.M1 = 0;
                this.M2 = 0;
                this.M3 = 0;
            end
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
            copy.M1 = this.M1;
            copy.M2 = this.M2;
            copy.M3 = this.M3;
            % No cloning of the associated reduced model.
            copy.ReducedModel = this.ReducedModel;
        end
        
        function set.ReducedModel(this, value)
            if ~isa(value,'models.ReducedModel')
                error('The given value has to be a models.ReducedModel instance.');
            end
            this.ReducedModel = value;
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
            demo = LocalLipEstimatorDemo(3);
            demo.Run;
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

