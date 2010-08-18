classdef ExperimentalEstimator < error.BaseEstimator
    %ExperimentalEstimator 
    %
    % Error estimator for the getCFI-attempt (non-gronwall-estimation)
    
    properties
        % Minimum threshold for `M_i` matrix values.
        %
        % All matrix entries from matrices `M_1,M_2,M_3` lower than
        % MinMatval are set to zero.
        % This is to allow reduction of numerical rounding errors.
        %
        % Default: eps
%        MinMatval = eps;
    end
    
    properties(Access=private)
        M1;
        M2;
        M3;
    end
    
    methods
        function this = ExperimentalEstimator(rmodel)
            
            % Validity checks
            msg = error.ExperimentalEstimator.validModelForEstimator(rmodel);
            if ~isempty(msg)
                error(msg);
            end
            
            % Call superclass constructor with model argument
            this = this@error.BaseEstimator(rmodel);
            this.ExtraODEDims = 1;
            
            % Only prepare matrices if projection was used
            if ~isempty(rmodel.V) && ~isempty(rmodel.W)
                % P = (I-VW^t)
                P = (eye(size(rmodel.V,1)) - rmodel.V*rmodel.W');
                D = P' * rmodel.G * P;
                
                % Get full d x N coeff matrix of approx
                Ma = rmodel.FullModel.Approx.Ma;
                
                hlp = Ma'*D*Ma;
                %hlp(hlp < this.MinMatval) = 0;
                this.M1 = hlp;
                
                % Only linear input conversion (B = const. matrix) allowed so
                % far! mu,0 is only to let
                if ~isempty(rmodel.FullModel.System.B)
                    try
                        B = rmodel.FullModel.System.B.evaluate([],[]);
                    catch ME
                        B = rmodel.FullModel.System.B.evaluate(0,rmodel.System.getRandomParam);
                        warning('Some:Id','Error estimator for current system will not work correctly! (B is not linear and mu-independent!');
                    end
                    % Filter too small entries
                    hlp = B'*D*B;
                    %hlp(hlp < this.MinMatval) = 0;
                    this.M2 = hlp;
                    hlp = Ma'*D*B;
                    %hlp(hlp < this.MinMatval) = 0;
                    this.M3 = hlp;
                end
                
                % STUFF WITH t,mu
                
                % Create a clone of the input conversion
                %Bp = rmodel.FullModel.System.B.project([],P);
                %this.M2 = @(t,mu) B(t,mu)'*D*B(t.mu);
                %this.M3 = Ma'*D*Bp;
                
            else
                % No projection means no error!
                this.M1 = 0;
                this.M2 = 0;
                this.M3 = 0;
            end
        end
        
        function e = evalODEPart(this, x, t, mu, ut)
            phi = this.ReducedModel.System.f.evaluateAtCenters(x(1:end-1), t, mu);            
            
            % I-WV^t part
            if nargin == 5
                e = phi'*this.M1*phi + ut'*this.M2*ut + phi'*this.M3*ut;
                % No input case
            else
                e = phi'*this.M1*phi;
            end
            e = sqrt(max(e,0));
            fci = this.ReducedModel.FullModel.System.f.getcfi(x(1:end-1), abs(x(end)), t, mu);
            e = e + fci;
        end
        
        function process(this, t, x, mu, inputidx)%#ok
            eint = x(end-this.ExtraODEDims+1:end,:);
            if all(eint == 0)
                warning('ExperimentalEstimator:process','Integral part is all zero. Attention!');
            end
            this.LastError = eint + this.ReducedModel.getExo(mu);
        end
        
        function e0 = getE0(this, mu)%#ok
            e0 = 0;
        end
        
    end
    
    methods(Static)
        function errmsg = validModelForEstimator(rmodel)
            % Validations
            errmsg = [];
            if ~isa(rmodel.System.f,'dscomponents.AKernelCoreFun')
                errmsg = 'The reduced model''s core function must be a subclass of dscomponents.AKernelCoreFun for this error estimator.'; 
            elseif ~isa(rmodel.FullModel.Approx,'dscomponents.CompwiseKernelCoreFun')
                errmsg = 'The full model''s approx function must be a subclass of dscomponents.CompwiseKernelCoreFun for this error estimator.'; 
            end
        end
    end
    
end