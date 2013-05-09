classdef ACoreFun < KerMorObject & general.AProjectable
% Basic interface for all dynamical system's core functions
% Inherits the AProjectable interface.
%
% Subclassers have to implement the clone method according to the rules
% explained in ICloneable, but only must re-implmement the project
% method if a custom projection process takes place, i.e. coefficients
% must be recomputed or such. Otherwise, the project implementation in
% this method is sufficient to wrap 
%
% @author Daniel Wirtz @date 2010-03-17
%
% @new{0,6,dw,2012-07-16} Removed the CustomJacobian flag as either the default via finite
% differences is used or an overridden method.
%
% @new{0,6,dw,2012-02-07}
% - Added a new property dscomponents.ACoreFun.CustomJacobian and
% implemented the default jacobian via finite differences. 
% - Removed the IJacobian interface and adopted the affected classes.
% - Implemented a test function to compare both jacobians (custom and FD)
%
% @change{0,6,dw,2012-01-18} Fixed a bug so that PointerCoreFuns with no `t` or `\mu` argument
% can be evaluated and no index out of range errors occur. Also the output dimension is now
% automatically determined in that case.
%
% @change{0,6,dw,2011-11-30} Fixed a bug regarding the setter for JSparsityPattern, as
% previously passing an empty value was not allowed.
%
% @new{0,6,dw,2011-11-27} New property ACoreFun.JSparsityPattern that allows to specify the
% sparsity pattern of the functions jacobian matrix. Is used for any solver that supports those
% patterns.
%
% @change{0,5,dw,2011-10-15} Improved the evaluate method and added a
% generic test method to test the multiargumentevaluation-capability (calls
% once with whole vector and then via singledim-loop, then compares diff up to eps)
%
% @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
% supervision system @ref propclasses. This class now inherits from KerMorObject and has an
% extended constructor registering any user-relevant properties using
% KerMorObject.registerProps.
%
% @change{0,3,sa,2011-04-15} Implemented Setters for the properties of this class 
%
% @change{0,3,dw,2011-04-13} The evaluate function supports multiargument evaluation now.
%
% @new{0,3,dw,2011-04-12} Added a new set-protected property
% dscomponents.ACoreFun.MultiArgumentEvaluations. Allows for speedup in case the target function
% can be called with a matrix of vectors instead of a single vector only.
%
% @todo Add test for non-custom projection (model with `W'f(x,t,\mu)V`)
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetObservable)
        % Flag that indicates if the ACoreFun is (truly) time-dependent.
        %
        % Set in subclasses to the correct value for your implementation.
        % Most built-in KerMor classes set their values correctly or guess them.
        %
        % This will cause different behaviour for e.g. ODE solvers.
        %
        % @propclass{critical} Not setting this value in implementing subclasses causes
        % KerMor's ODE solvers to (possibly) produce wrong results due to wrong assumptions on
        % the time dependence of the core function.
        %
        % @type logical @default true
        TimeDependent = true;
    end
    
    properties(SetObservable, SetAccess=protected)
        % Set this property if the projection process is customized by
        % overriding the default project method.
        % 
        % This affects the evaluation method of the ACoreFun.
        %
        % @propclass{optional} If the structure of the underlying function allows more efficient
        % projection set this flag to true in the subclass constructor.
        CustomProjection = false;
        
        % Flag that indicates whether this core function can be evaluated passing a matrix of
        % vectors instead of a single vector
        %
        % In case the custom implementation allows for such a call a significant speedup can be
        % achieved when using matlabs vector-based expressions
        %
        % @note Regarding parameters and their processing in multiargument-evaluations, the bsxfun
        % method from matlab may be of great help. See the models.pcd.CoreFun1D class for an
        % example.
        %
        % @propclass{optional} Speedup can be gained if subclasses allow matrices to be passed
        % instead of single vectors.
        %
        % @default false
        %
        % See also: evaluate bsxfun models.pcd.CoreFun1D
        MultiArgumentEvaluations = false;
        
        % Sparsity pattern for the jacobian matrix.
        %
        % KerMor automatically detects nonempty values and forwards them as appropriate to any
        % ODE solver.
        %
        % @propclass{optional} Some ODE solvers can work more efficiently if a sparsity pattern
        % for the jacobian matrix of the core function can be provided.
        % 
        % @type sparse<logical> @default []
        %
        % @todo maybe move this property to the IJacobian interface? (but:
        % might have sparsity pattern but not actually a analytic
        % jacobian..)
        JSparsityPattern = [];
        
        % The current state space dimension of the function's argument `x`.
        %
        % Must be set in inheriting classes on order to ensure
        % functionality of KerMor classes and algorithms.
        %
        % @propclass{critical} Not setting this value in inheriting classes
        % may cause errors in some KerMor algorithms.
        %
        % @type integer @default []
        xDim = [];
        
        % The current output dimension of the function.
        %
        % Must be set in inheriting classes on order to ensure
        % functionality of KerMor classes and algorithms.
        %
        % @propclass{critical} Not setting this value in inheriting classes
        % may cause errors in some KerMor algorithms.
        %
        % @type integer @default []
        fDim = [];
    end
    
    methods
        
        function this = ACoreFun
            this = this@KerMorObject;
            this.registerProps('CustomProjection',...
                'MultiArgumentEvaluations','JSparsityPattern',...
                'TimeDependent','xDim','fDim');
        end
        
        function target = project(this, V, W, target)
            % Sets the protected `\vV,\vW` matrices that can be utilized on core
            % function evaluations after projection.
            %
            % If the customized subclass sets the CustomProjection property
            % to true, this method must be called with a 3rd argument, the
            % target instance of a subclass upon which the projection
            % process is to be performed.
            %
            % Parameters:
            % V: Projection matrix one
            % W: Projection matrix two
            % target: The target instance which is a clone of this current
            % class. This parameter is required if custom projection is
            % used.
            %
            % See also: CustomProjection
            %
            % @change{0,3,dw,2011-04-11} Modified the project method to
            % also cater for the noncustomized projection case by cloning
            % the current instance (general.AProjectable now inherits
            % from ICloneable per default).
            if nargin < 4
                if this.CustomProjection
                    error('The target parameter must be given if custom projection is used.');
                end
                target = this.clone;
            end
            target = project@general.AProjectable(this, V, W, target);
            % New state space argument size is the number of columns of V.
            if isempty(V)
                target.xDim = this.xDim;
            else
                target.xDim = size(V,2);
            end
            % New state space output size is the number of columns of W.
            if isempty(W)
                target.fDim = this.fDim;
            else
                target.fDim = size(W,2);
            end
        end
        
        function fx = evaluate(this, x, t, mu)
            % Evaluates the f-approximation. Depending on a possible
            % projection and the CustomProjection-property the function
            % either calls the inner evaluation directly which assumes
            % `f = f^r(z)` or projects the reduced state variable z into
            % the original space and evaluates the function there, so via
            % `f = V'f(Vz)`
            %
            % Attention:
            % If you override this method directly, you will have to make sure it can handle
            % matrix-based `x` inputs and causes correct evaluation after projection has been
            % performed.
            %
            % Parameters:
            % x: The state variable vector/matrix (with colum state
            % vectors)
            % t: The corresponding times for each state vector. Set to []
            % if no time is used.
            % mu: The parameter(s) to use. Set to [] if the function does not
            % support parameters.
            %
            % See also: CustomProjection
            %
            % @change{0,3,dw,2011-04-19} Fixed an error when no MultiArgumentEvaluations were supported and
            % no parameter `\mu` was given (Crashed due to index out of bounds)
            
            proj = ~this.CustomProjection && ~isempty(this.V) && ~isempty(this.W);
            if proj
                x = this.V*x;
            end
            % check if fast evaluation is possible
            if size(x,2) == 1 || this.MultiArgumentEvaluations
                fx = this.evaluateCoreFun(x, t, mu);
            else
                % evaluate each point extra
                if isempty(mu)
                    mu = double.empty(0,size(x,2));
                end
                if isempty(t)
                    t = double.empty(0,size(x,2));
                end
                % Detect output size from first evaluation
                fx = zeros(this.fDim, size(x,2));
                for idx = 1:size(x,2)
                    fx(:,idx) = this.evaluateCoreFun(x(:,idx), t(:,idx), mu(:,idx));
                end    
            end
            if proj
                fx = this.W'*fx;
            end
        end

        function J = getStateJacobian(this, x, t, mu)
            % Default implementation of jacobian matrix evaluation via
            % finite differences.
            %
            % Override in subclasses for custom (analytic) implementations or if the problem
            % size is too large to fit a full matrix into memory.
            %
            % Parameters:
            % x: The state variable vector/matrix (with colum state
            % vectors)
            % t: The corresponding times for each state vector. Set to []
            % if no time is used.
            % mu: The parameter(s) to use. Set to [] if the function does not
            % support parameters.
            %
            % Return values:
            % J: The jacobian matrix at `x,t,\mu`.
            dt = sqrt(eps);
            d = size(x,1);
            
            X = repmat(x,1,d); T = repmat(t,1,d); MU = repmat(mu,1,d);
            I = speye(d,d)*dt;
            % Evaluate makes use of built-in multi-argument evaluation
            J = (this.evaluate(X+I,T,MU) - this.evaluate(X,T,MU))/dt;
            % Create sparse matrix if pattern is set
            if ~isempty(this.JSparsityPattern)
                [i, j] = find(this.JSparsityPattern);
                J = sparse(i,j,J(logical(this.JSparsityPattern)),this.fDim,this.xDim);
            end            
        end
        
        function copy = clone(this, copy)
            if nargin == 1 || ~isa(copy,'dscomponents.ACoreFun')
                error('Incorrect call to clone. As this class is abstract, a subclass of ACoreFun has to be passed as second argument.');
            end
            copy = clone@general.AProjectable(this, copy);
            % Copy local properties
            copy.CustomProjection = this.CustomProjection;
            copy.MultiArgumentEvaluations = this.MultiArgumentEvaluations;
            copy.JSparsityPattern = this.JSparsityPattern;
            copy.TimeDependent = this.TimeDependent;
            copy.xDim = this.xDim;
            copy.fDim = this.fDim;
        end
    end
        
    methods(Abstract)
        % Actual method used to evaluate the dynamical sytems' core function.
        %
        % Subclasses might implement this method and set the flags CustomProjection and
        % MultiArgumentEvaluations appropriately.
        % However, for speed reasons, if both are true one might as well override the 'evaluate' member directly as it
        % basically cares for the cases when one of the flags is not true. In that case it is still
        % important to set both flags to true as some components rely on them.
        fx = evaluateCoreFun(this, x, t, mu);
    end
    
    %% Getter & Setter
    methods
        function set.CustomProjection(this, value)
            if ~islogical(value)
                error('Property must be logical/boolean. either true or false');
            end
            this.CustomProjection = value;
        end
        
        function set.MultiArgumentEvaluations(this, value)
            if ~islogical(value)
                error('Property must be logical/boolean. either true or false');
            end
            this.MultiArgumentEvaluations = value;
        end
        
        function set.JSparsityPattern(this, value)
            if ~isempty(value) && ~issparse(value)
                error('JSparsityPattern must be a sparse matrix.');
            end
            this.JSparsityPattern = value;
        end
        
        function set.xDim(this, value)
            if ~isempty(value) && ~isposintscalar(value)
                error('xDim must be a positive integer.');
            end
            this.xDim = value;
        end
        
        function set.fDim(this, value)
            if ~isempty(value) && ~isposintscalar(value)
                error('fDim must be a positive integer.');
            end
            this.fDim = value;
        end
        
        function res = test_MultiArgEval(this, mudim)
            % Convenience function that tests if a custom
            % MultiArgumentEvaluation works as if called with single
            % arguments.
            %
            if nargin < 2
                mudim = 100;
            end
            if true || (this.MultiArgumentEvaluations)
                x = rand(this.xDim,200);
                mu = rand(mudim,200);
                fxm = this.evaluate(x,1:200,mu);
                fxs = zeros(size(x));
                for i=1:200
                    fxs(:,i) = this.evaluate(x(:,i),i,mu(:,i));
                end
                err = sum((fxm-fxs).^2,1);
                res = all(err < eps);
                if ~res
                    plot(err);
                end
            else
                error('MultiArgumentEvaluations not switched on.');
            end
        end
        
        function res = test_Jacobian(this, xa, ta, mua)
            % Tests the custom provided jacobian matrix against the default
            % finite difference computed one.
            %
            % Attention:
            % This method trivially works if you did not overload the
            % dscomponents.ACoreFun.getStateJacobian method (then it will check the default
            % implementation against itself)
            %
            % Parameters:
            % xa: A matrix of `x_i,i=1\ldots n`-values to try each. @type matrix<double>
            % ta: The corresponding `n` times `t_i` @type rowvec<double>
            % mua: The parameters `\mu_i`. @type matrix<double> @default []
            %
            % Return values:
            % res: A flag indicating if each test had a relative error of
            % less than 1e-7 @type logical
            reltol = 1e-7;
            dt = sqrt(eps);
            res = false;
            if nargin < 3 || isempty(ta)
                ta = double.empty(0,size(xa,2));
            end
            if nargin < 4 || isempty(mua)
                mua = double.empty(0,size(xa,2));
            end
            d = size(xa,1);
            I = speye(d,d)*dt;
            for i = 1:size(xa,2)
                x = xa(:,i); t = ta(:,i); mu = mua(:,i);
                Jc = this.getStateJacobian(x, t, mu);
                
                %% Numerical jacobian
                perstep = floor((256*1024^2)/(8*d));
                X = repmat(x,1,perstep); T = repmat(t,1,perstep); MU = repmat(mu,1,perstep);
                steps = ceil(d/perstep);
                pi = ProcessIndicator('Comparing %dx%d jacobian with finite differences over %d blocks of size %d',...
                    steps,false,this.fDim,this.xDim,steps,perstep);
                for k = 1:steps
                    if k == steps
                        num = d-(k-1)*perstep;
                        X = repmat(x,1,num); T = repmat(t,1,num); MU = repmat(mu,1,num);
                    end
                    pos = (k-1)*perstep+1:min(d,k*perstep);
                    J = (this.evaluate(X+I(:,pos),T,MU) - this.evaluate(X,T,MU))/dt;
                    pi.step;
                    abserr = max(max(abs(J-Jc(:,pos))));
                    relerr = abserr/max(max(abs(J)));
                    if relerr >= reltol;
                        pi.stop;
                        fprintf('Failed. Max absolute error %g, relative %g\n',abserr,relerr);
                        return;
                    end
                end
                pi.stop;
            end
            res = true;
        end
    end
end

