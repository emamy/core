classdef BaseLipKernelEstimator < error.BaseEstimator
    %BASELIPKERNELESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private, GetAccess=protected)
        M1 = [];
        M2 = [];
        M3 = [];
    end
    
    methods
        function setReducedModel(this, rmodel)
            % Overrides the method from BaseEstimator and performs
            % additional computations.
            
            % Call superclass
            setReducedModel@error.BaseEstimator(this, rmodel);
            
            % Perform any offline computations/preparations
            % Only prepare matrices if projection is used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                
                fm = rmodel.FullModel;
                
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
                
                % Compute projection part matrices, without creating a
                % d x d matrix (too big!)
                M = rmodel.V*(rmodel.W'*Ma);
                G1 = Ma'*rmodel.GScaled;
                this.M1 = G1*Ma - 2*G1*M + M'*(rmodel.GScaled*M);
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(fm.System.B)
                    try
                        B = fm.System.B.evaluate([],[]);
                    catch ME%#ok
                        B = fm.System.B.evaluate(0,rmodel.System.getRandomParam);
                        warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                    end
                    
                    B2 = rmodel.V*(rmodel.W'*B);
                    G2 = B'*rmodel.GScaled;
                    this.M2 = 2*(G1*B - M'*G2' - G1*B2 + M'*(rmodel.GScaled*B2));
                    this.M3 = G2*B - 2*G2*B2 + B2'*(rmodel.GScaled*B2);
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
        
        function copy = clone(this, copy)
            copy = clone@error.BaseEstimator(this, copy);
            copy.M1 = this.M1;
            copy.M2 = this.M2;
            copy.M3 = this.M3;
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
            end
        end
    end
    
end

