classdef BaseCompLemmaEstimator < error.BaseEstimator
    % BaseCompLemmaEstimator: Base class for error estimators using the comparison lemma formulation.
    %
    % @author Daniel Wirtz @date 2010-08-10
    %
    % @change{0,5,dw,2011-07-04} Changed this classes name to "BaseCompLemmaEstimator".
    %
    % @change{0,4,dw,2011-05-29} Changed this classes name to "BaseKernelEstimator".
    %
    % @new{0,1,dw,2010-08-10} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    %
    % @todo make listener fields transient and re-register upon loading!
    %
    % @todo remove old ISimConstants interface implementation
    % prepareConstants, at least rename it somehow
    
    properties(Transient, SetAccess=private, GetAccess=protected)
        % The current step
        StepNr;
    end
    
    properties(SetAccess=protected)
        % `3\times n` matrix containing the `t,\alpha(t),\beta(t)` values at each time step.
        EstimationData = [];
    end
    
    properties(SetAccess=private)
        aComp;
    end
    
    properties(Access=private, Transient)
        % The cbPreSolve listener instance
        lstPreSolve;
    end
        
    methods
        function this = BaseCompLemmaEstimator
            this = this@error.BaseEstimator;
            this.ExtraODEDims = 1;
        end
        
        function offlineComputations(this, model)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            
            % Call superclass
            offlineComputations@error.BaseEstimator(this, model);
            
            if isa(model.System.B,'dscomponents.AffLinInputConv')
                this.aComp = error.alpha.AffineParametric(model);
            else
                this.aComp = error.alpha.Constant(model);
            end
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            % Evaluates the auxiliary ode part for the comparison-lemma
            % based error estimators.
            %
            % Parameters:
            % x: The full extended state variable vector. Extended means that
            % the last @ref ExtraODEDims rows contain the error estimators own
            % data. If not used, implementers must take care to ditch those
            % values if any function evaluations are performed within the
            % integral part. @type colvec
            % t: The current time `t` @type double
            % mu: The current parameter `\mu` @type colvec
            % ut: The value of the input function `u(t)` if given, [] else.
            % @type double
            %
            % Return values:
            % e: The auxiliary ode part value.
            phi = this.ReducedModel.System.f.getKernelVector(x(1:end-1), t, mu);
            
            a = this.aComp.getAlpha(phi, ut, t, mu);
            b = this.getBeta(x, t, mu);
            e = b*x(end) + a;
            
            this.EstimationData(:,this.StepNr) = [t; a; b];
            this.StepNr = this.StepNr + 1;
        end
        
        function ct = prepareConstants(this, ~, ~)
            st = tic;
            if isempty(this.lstPreSolve)
                this.lstPreSolve = addlistener(rmodel.ODESolver,'PreSolve',@this.cbPreSolve);
            end
            this.lstPreSolve.Enabled = true;
            this.StepNr = 1;
            ct = toc(st);
        end
        
        function copy = clone(this, copy)
            copy = clone@error.BaseEstimator(this, copy);
            copy.EstimationData = this.EstimationData;
            copy.aComp = this.aComp;
            copy.StepNr = this.StepNr;
            copy.lstPreSolve = addlistener(copy.ReducedModel.ODESolver,'PreSolve',@copy.cbPreSolve);
            copy.lstPreSolve.Enabled = this.lstPreSolve.Enabled;
        end
        
        function clear(this)
            clear@error.BaseEstimator(this);
            this.EstimationData = [];
        end
    end
    
    methods(Access=protected)
        function ct = postprocess(this, x, t, mu, inputidx)%#ok
            % Return values:
            % ct: The time needed for postprocessing @type double
            if ~isempty(this.lstPreSolve)
                this.lstPreSolve.Enabled = false;
            end
            ct = 0;
        end
    end
    
    methods(Abstract, Access=protected)
        % Computes the `\beta(t)` term from the error estimation ODE for given time and place.
        %
        % Parameters:
        % x: The current reduced state variable, composed by `Vz(t)` and any extra dimensions set up
        % by the error estimator. @type colvec
        % t: The current time `t\in[0,T]` @type double
        % mu: The current parameter `\mu` @type colvec
        b = getBeta(this, x, t, mu);
    end
    
    methods(Access=private)
        function cbPreSolve(this, sender, data)%#ok
            % Ensure that the estimation data matrix has the correct size but leave it as-is for
            % possibly having iterations in subclasses.
            if isempty(this.EstimationData) || size(this.EstimationData,2) ~= length(data.Times)
                this.EstimationData = zeros(3,length(data.Times));
            end
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(model)
            % Validations
            errmsg = [];
            if ~isempty(model.Approx) 
                if ~isa(model.Approx,'approx.KernelApprox')
                    errmsg = 'The full model''s approx function must be a subclass of approx.KernelApprox for this error estimator.'; 
                    return;
                end
            elseif ~isa(model.System.f,'kernels.ParamTimeKernelExpansion')
                    errmsg = 'If no approximation is used, the full model''s core function must be a subclass of kernels.ParamTimeKernelExpansion for this error estimator.'; 
                    return;
            end
            if ~isa(model.System.C,'dscomponents.LinearOutputConv')
                errmsg = 'Local Lipschitz estimators work only for constant linear output conversion.';
                return;
            elseif model.System.C.TimeDependent
                errmsg = 'Output error estimation for time dependent output not implemented yet.';
                return;
            end
            if ~isa(model.ODESolver,'solvers.BaseCustomSolver');
                errmsg = 'The reduced models ODE solver must be a subclass of BaseCustomSolver.';
            end
        end
    end
end

