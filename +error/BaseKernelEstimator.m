classdef BaseKernelEstimator < error.BaseEstimator
    % BaseKernelEstimator: Base class for local lipschitz error estimators.
    %
    % @author Daniel Wirtz @date 2010-08-10
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
        
    properties(SetAccess=private, GetAccess=protected)
        M1 = [];
        M2 = [];
        M3 = [];
    end
    
    properties(Transient, SetAccess=private, GetAccess=protected)
        % The current step
        StepNr;
    end
    
    properties(SetAccess=protected)
        % `3\times n` matrix containing the `t,\alpha(t),\beta(t)` values at each time step.
        EstimationData = [];
    end
    
    properties(Access=private, Transient)
        % The cbPreSolve listener instance
        lstPreSolve;
    end
        
    methods
        function this = BaseKernelEstimator
            this = this@error.BaseEstimator;
            this.ExtraODEDims = 1;
        end
        
        function setReducedModel(this, rmodel)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            
            % Call superclass
            setReducedModel@error.BaseEstimator(this, rmodel);
            
            fm = rmodel.FullModel; 
            B = [];
            if ~isempty(fm.System.B)
                try
                    B = fm.System.B.evaluate([],[]);
                catch ME%#ok
                    B = fm.System.B.evaluate(0,rmodel.System.getRandomParam);
                    warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                end
            end
            
            % Obtain the correct snapshots
            % Standard case: the approx function is a kernel expansion. it
            % can also be that the system's core function is already a
            % kernel expansion
            if ~isempty(fm.Approx)
                % Get full d x N coeff matrix of approx function
                Ma = fm.Approx.Ma;
            else
                % Get full d x N coeff matrix of core function
                Ma = fm.System.f.Ma;
            end
            
            % Perform any offline computations/preparations
            % Only prepare matrices if projection is used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                
                % Compute projection part matrices, without creating a
                % d x d matrix (too big!)
                M = Ma - rmodel.V*(rmodel.W'*Ma);
                hlp = M'*(rmodel.GScaled*M);
                % Check if matrix needs to be made symmetric
                if any(any(abs(hlp-hlp') > 1e-5))
                    hlp = (hlp + hlp')/2;
                    warning('KerMor:errorest','M1 matrix not sufficiently symmetric, updating (M+M'')/2');
                end
                this.M1 = hlp;
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(B)
                    B2 = B-rmodel.V*(rmodel.W'*B);
                    this.M2 = M'*(rmodel.GScaled*B2);
                    this.M3 = B2'*(rmodel.GScaled*B2);
                    clear B2;
                end
                clear M;
            else
                % No projection means no projection error!
                n = size(Ma,2);
                this.M1 = zeros(n,n);
                if ~isempty(B)
                    b = size(B,2);
                    this.M2 = zeros(n,b);
                    this.M3 = zeros(b,b);
                end
            end
            
            this.lstPreSolve = addlistener(rmodel.ODESolver,'PreSolve',@this.cbPreSolve);
            this.lstPreSolve.Enabled = false;
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            
            % Compute `\alpha(t)`
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-1), t, mu);
            
            a = phi*this.M1*phi';
            if ~isempty(ut) % An input function u is set
                a = a + phi*this.M2*ut + ut'*this.M3*ut;    
            end
            a = sqrt(abs(a));
            %a = sqrt(max(a,0));
            
            b = this.getBeta(x, t, mu);
            e = b*x(end) + a;
            
            this.EstimationData(:,this.StepNr) = [t; a; b];
            this.StepNr = this.StepNr + 1;
        end
        
        function prepareConstants(this)
            this.lstPreSolve.Enabled = true;
            this.StepNr = 1;
        end
        
        function e0 = getE0(this, mu)
            % Returns the initial error at `t=0` of the integral part.
            e0 = this.ReducedModel.getExo(mu);
        end
        
        function copy = clone(this, copy)
            copy = clone@error.BaseEstimator(this, copy);
            copy.M1 = this.M1;
            copy.M2 = this.M2;
            copy.M3 = this.M3;
            copy.EstimationData = this.EstimationData;
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
        function postprocess(this, t, x, mu, inputidx)%#ok
            this.lstPreSolve.Enabled = false;
        end
    end
    
    methods(Abstract, Access=protected)
        % Computes the `\beta(t)` term from the error estimation ODE for given time and place.
        %
        % Parameters:
        % x: The current reduced state `Vz(t)`
        % t: The current time `t\in[0,T]`
        % mu: The current parameter `\mu`
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
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = [];
            if ~isempty(rmodel.FullModel.Approx) 
                if ~isa(rmodel.FullModel.Approx,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'The full model''s approx function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
                    return;
                end
            elseif ~isa(rmodel.FullModel.System.f,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'If no approximation is used, the full model''s core function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
                    return;
            end
            if ~isa(rmodel.FullModel.System.C,'dscomponents.LinearOutputConv')
                errmsg = 'Local Lipschitz estimators work only for constant linear output conversion.';
                return;
            elseif rmodel.FullModel.System.C.TimeDependent
                errmsg = 'Output error estimation for time dependent output not implemented yet.';
                return;
            end
            if ~isa(rmodel.ODESolver,'solvers.ode.BaseCustomSolver');
                errmsg = 'The reduced models ODE solver must be a subclass of BaseCustomSolver.';
            end
        end
    end
    
     methods(Static,Access=protected)
        function s = loadobj(s)
            s = loadobj@KerMorObject(s);
            s.lstPreSolve = addlistener(s.ReducedModel.ODESolver,'cbPreSolve',@s.cbPreSolve);
            s.lstPreSolve.Enabled = false;
        end
    end
    
end

