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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        % @todo maybe move this property to the IJacobian interface?
        % (but: might have sparsity pattern but not actually a analytic
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
    
    properties(SetAccess=private)
        % The system associated with the current ACoreFun
        %
        % @type models.BaseDynSystem
        System;
    end
    
    properties(SetAccess=private, Transient)
        % The current model parameter mu for evaluations. Will not be
        % persisted as only valid for runtime during simulations.
        %
        % @type colvec<double> @default []
        %
        % See also evaluate
        mu = [];
    end
    
    methods
        
        function this = ACoreFun(sys)
            this = this@KerMorObject;
            this.registerProps('CustomProjection',...
                'MultiArgumentEvaluations','JSparsityPattern',...
                'TimeDependent','xDim','fDim');
            this.System = sys;
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
        
        function fx = evaluate(this, x, t)
            % Evaluates the f-approximation. Depending on a possible
            % projection and the CustomProjection-property the function
            % either calls the inner evaluation directly which assumes
            % `f = f^r(z)` or projects the reduced state variable z into
            % the original space and evaluates the function there, so via
            % `f = V'f(Vz)`
            %
            % Attention:
            % If you override this method directly, you will have to make
            % sure it computes the correct evaluation after projection has
            % been performed.
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
            fx = this.evaluateCoreFun(x, t);
            if proj
                fx = this.W'*fx;
            end
        end
        
        function fx = evaluateMulti(this, x, t, mu)
            % Evaluates this function on multiple locations and maybe
            % multiple times and parameters.
            %
            % The parameters t and mu can either be single or as many as
            % there are columns in x. Overriding classes must adhere to
            % this principle.
            %
            % Override in subclasses if the "evaluate" can handle matrix
            % arguments directly
            singlemu = size(mu,2) == 1;
            if isempty(mu)
                mu = [];
                singlemu = true;
            end
            if singlemu
                this.prepareSimulation(mu);
            end
            if isempty(t)
                t = double.empty(0,size(x,2));
            elseif length(t) == 1
                t = ones(1,size(x,2))*t;
            end
            m = size(x,2);
            if m == 0
                m = length(t);
            end
            fx = zeros(this.fDim, m);
            for idx = 1:size(x,2)
                if ~singlemu
                    this.prepareSimulation(mu(:, idx));
                end
                fx(:,idx) = this.evaluate(x(:,idx), t(:,idx));
            end
        end
        
        function prepareSimulation(this, mu)
            % A method that allows parameter-dependent computations to be
            % performed before a simulation using this parameter starts.
            %
            % Override in subclasses (with calling this method!) for
            % customized behaviour
            this.mu = mu;
        end

        function J = getStateJacobian(this, x, t)
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
            J = this.getStateJacobianFD(x, t);
        end
        
        function copy = clone(this, copy)
            if nargin == 1 || ~isa(copy,'dscomponents.ACoreFun')
                error('Incorrect call to clone. As this class is abstract, a subclass of ACoreFun has to be passed as second argument.');
            end
            copy = clone@general.AProjectable(this, copy);
            % Copy local properties
            copy.CustomProjection = this.CustomProjection;
            % Dont clone the associated system
            copy.System = this.System;
            copy.JSparsityPattern = this.JSparsityPattern;
            copy.TimeDependent = this.TimeDependent;
            copy.xDim = this.xDim;
            copy.fDim = this.fDim;
        end
    end
    
    methods(Sealed)
        function J = getStateJacobianFD(this, x, t, partidx)
            % Implementation of jacobian matrix evaluation via
            % finite differences.
            %
            % Parameters:
            % x: The state variable vector/matrix (with colum state
            % vectors)
            % t: The corresponding times for each state vector. Set to []
            % if no time is used.
            % mu: The parameter(s) to use. Set to [] if the function does not
            % support parameters.
            % partidx: An index vector for desired derivative components.
            % Can be used to compute full jacobians in parts that would
            % otherwise not fit into memory @type rowvec<integer> @default
            % 1:xDim
            %
            % Return values:
            % J: The jacobian matrix at `x,t,\mu`, possibly only the part
            % specified by the partidx indices. @type matrix<double>
            dt = sqrt(eps);
            d = size(x,1);
            if nargin < 4
                partidx = 1:this.xDim;
            end
            len = length(partidx);
            X = repmat(x,1,len);
            I = speye(d,d)*dt;
            del = 1:this.xDim;
            del(partidx) = [];
            I(:,del) = [];
            % Evaluate makes use of built-in multi-argument evaluation
            J = (this.evaluateMulti(X+I,t,this.mu) - repmat(this.evaluate(x,t),1,len))/dt;
            % Create sparse matrix if pattern is set
            if ~isempty(this.JSparsityPattern) && nargin < 4
                [i, j] = find(this.JSparsityPattern);
                J = sparse(i,j,J(logical(this.JSparsityPattern)),this.fDim,this.xDim);
            end            
        end
    end
        
    methods(Abstract)
        % Actual method used to evaluate the dynamical sytems' core function.
        %
        % Subclasses might implement this method and set the flag
        % CustomProjection appropriately.
        % However, for speed reasons, if both are true one might as well
        % override the 'evaluate' member directly as it
        % basically cares for the cases when one of the flags is not true.
        % In that case it is still
        % important to set both flags to true as some components rely on
        % them.
        fx = evaluateCoreFun(this, x, t);
    end
    
    %% Getter & Setter
    methods
        function set.CustomProjection(this, value)
            if ~islogical(value)
                error('Property must be logical/boolean. either true or false');
            end
            this.CustomProjection = value;
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
            x = rand(this.xDim,200);
            mui = rand(mudim,200);
            fxm = this.evaluateMulti(x,1:200,mui);
            fxs = zeros(size(x));
            for idx = 1:size(x,2)
                this.prepareSimulation(mui(:, idx));
                fxs(:,idx) = this.evaluate(x(:,idx), idx);
            end
            err = sum((fxm-fxs).^2,1);
            res = all(err < eps);
            if ~res
                plot(err);
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
            if nargin < 2
                xa = rand(this.xDim,20);
                ta = 1:20;
                mua = rand(20,20);
            else
                if nargin < 4 || isempty(mua)
                    mua = double.empty(0,size(xa,2));
                    if nargin < 3 || isempty(ta)
                        ta = double.empty(0,size(xa,2));
                    end
                end
            end
            d = size(xa,1);
            perstep = floor((256*1024^2)/(8*d)); % 256MB chunks
            steps = ceil(d/perstep);
            gpi = steps == 1 && size(xa,2) > 1;
            if gpi
                pi = ProcessIndicator('Comparing %d %dx%d jacobians with finite differences',...
                    size(xa,2),false,size(xa,2),this.fDim,this.xDim);
            end
            if size(mua,2) == 1
                mua = repmat(mua,1,size(xa,2));
            end
            for i = 1:size(xa,2)
                x = xa(:,i); t = ta(:,i);
                this.prepareSimulation(mua(:,i))
                Jc = this.getStateJacobian(x, t);
                
                %% Numerical jacobian
                if ~gpi
                    pi = ProcessIndicator('Comparing %dx%d jacobian with finite differences over %d blocks of size %d',...
                        steps,false,this.fDim,this.xDim,steps,perstep);
                end
                for k = 1:steps
                    num = 6;
%                     if k == steps
%                         num = d-(k-1)*perstep;
%                     end
                    pos = (k-1)*perstep+1:min(d,k*perstep);
                    J = this.getStateJacobianFD(x, t, pos);
                    pi.step;
                    abserr = max(max(abs(J-Jc(:,pos))));
                    relerr = max(max(abs(J-Jc(:,pos))./J));
                    [posi,posj] = ind2sub(size(J), pos(1:num));
                    indic = sprintf('%d, ',posi(1:num));
                    indjc = sprintf('%d, ',posj(1:num));
                    if relerr >= reltol;
                        diff = abs(J-Jc(:,pos));
                        [v, pos] = sort(diff(:),'descend');
                        
                        maxJ = max(max(abs(J)));
                        M = [v(1:num) v(1:num)./J(pos(1:num)) v(1:num)/maxJ];
                        
                        pi.stop;
                        disp(M);
                        fprintf('Failed at test vector %d. Max absolute error %g, relative %g, max %d errors at rows %s, cols %s (maxJ=%g)\n',i,abserr,relerr,num,indic,indjc,J(pos(1)));
                        %return;
                    else
                        fprintf('Max absolute error %g, relative %g, max %d errors at rows %s, cols %s (maxJ=%g)\n',abserr,relerr,num,indic,indjc,J(pos(1)));
                    end
                end
                if ~gpi
                    pi.stop;
                end
            end
            pi.stop;
            res = true;
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if nargin == 2
                obj.CustomProjection = from.CustomProjection;
                obj.JSparsityPattern = from.JSparsityPattern;
                obj.TimeDependent = from.TimeDependent;
                obj.xDim = from.xDim;
                obj.fDim = from.fDim;
                obj = loadobj@general.AProjectable(obj, from);
                obj = loadobj@KerMorObject(obj, from);
            else
                obj = loadobj@general.AProjectable(obj);
                obj = loadobj@KerMorObject(obj);
            end
        end
    end
end

