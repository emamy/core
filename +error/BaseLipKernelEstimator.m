classdef BaseLipKernelEstimator < error.BaseEstimator
    % BaseLipKernelEstimator: Base class for local lipschitz error estimators.
    %
    % @author Daniel Wirtz @date 2010-08-10
    %
    % @new{0,1,dw,2010-08-10} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(SetAccess=private, GetAccess=protected)
        M1 = [];
        M2 = [];
        M3 = [];
    end
    
    properties(SetAccess=protected)
        betas = [];
    end
    
    methods
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
        end
        
        function copy = clone(this, copy)
            copy = clone@error.BaseEstimator(this, copy);
            copy.M1 = this.M1;
            copy.M2 = this.M2;
            copy.M3 = this.M3;
        end
        
        function clear(this)
            clear@error.BaseEstimator(this);
            this.betas = [];
        end
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = [];
            if ~isempty(rmodel.FullModel.Approx) 
                if ~isa(rmodel.FullModel.Approx,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'The full model''s approx function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
                end
            elseif ~isa(rmodel.FullModel.System.f,'dscomponents.CompwiseKernelCoreFun')
                    errmsg = 'If no approximation is used, the full model''s core function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
            end
            if ~isa(rmodel.FullModel.System.C,'dscomponents.LinearOutputConv')
                errmsg = 'Local Lipschitz estimators work only for constant linear output conversion.';
            elseif rmodel.FullModel.System.C.TimeDependent
                errmsg = 'Output error estimation for time dependent output not implemented yet.';
            end
        end
    end
    
end

