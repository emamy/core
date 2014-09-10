classdef ACompEvalCoreFun < dscomponents.ACoreFun
% ACompEvalCoreFun: A normal CoreFun which supports single-component
% evaluation.
%
% This extra abilities are used within the approx.DEIM approximation.
%
% @author Daniel Wirtz @date 2012-03-26
%
% @change{0,6,dw,2012-06-11} The evaluateComponentPartialDerivatives method
% now supports vectorized inputs.
%
% @new{0,6,dw,2012-06-06} Added new methods
% ACompEvalCoreFun.evaluateJacobianSet and
% ACompEvalCoreFun.evaluateComponentPartialDerivatives to support jacobian
% entry evaluation for this class of functions (they inherit the same
% sparsity access to components). Also updated the test_CompEvalMatch to
% use random point sets and not only the full set (but that also)
%
% @new{0,6,dw,2012-05-30} Added the
% ACompEvalCoreFun.evaluateComponentSetGradientsAt and
% ACompEvalCoreFun.evaluateComponentGradientsAt methods for efficient
% jacobian computation of the DEIM approximation. The
% evaluateComponentGradientsAt default implementation uses finite
% differences.
%
% @new{0,6,dw,2012-03-26} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Dependent)
        PointSets;
    end
    
    properties(Access=private)
        pts = {};
        
        % Notational vectors for EI-point function argument dependencies
        jrow = {};
        jend = {};
        jself = {};
        deriv = {};
        dfxsel = {};
        T = {};
    end
    
    properties(SetAccess=private, GetAccess=protected)
        % The x-component selection matrices (precomputed on setting
        % PointSet/AltPointSet). �S\vx� is passed to the function
        % evaluating the components of �\vf� or its derivatives.
        %
        % See also: PointSet
        %
        % cellarray of @type matrix<double>
        S = {};
    end

    methods
        
        function this = ACompEvalCoreFun(sys)
            this = this@dscomponents.ACoreFun(sys);
        end
        
        function fx = evaluateComponentSet(this, nr, x, t)
            % Computes the full or reduced component functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            %
            % See also: PointSet
            fx = this.evaluateComponents(this.PointSets{nr},...
                this.jend{nr}, this.jrow{nr}, this.jself{nr}, this.S{nr}*x, t);
        end

        function fx = evaluateComponentSetMulti(this, nr, x, t, mu)
            % Computes the full component functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space locations `\vx` @type matrix<double>
            % t: The corresponding times `t` for the state `\vx` @type
            % rowvec<double>
            % mu: The corresponding parameters `\mu`. Can be a single
            % parameter or as many as the size of x @type matrix<double>
            %
            % See also: PointSet
            fx = this.evaluateComponentsMulti(this.PointSets{nr},...
                this.jend{nr}, this.jrow{nr}, this.jself{nr}, this.S{nr}*x, t, mu);
        end
        
        function dfx = evaluateComponentSetGradientsAt(this, nr, x, t)
            % Computes the full/reduced gradients of all component
            % functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            %
            % See also: PointSet
            %
            % Uses the template method, whose default implementation is via
            % finite differences.
            dfx = this.evaluateComponentGradientsAt(this.PointSets{nr},...
                    this.jend{nr}, this.jrow{nr}, this.jself{nr},...
                    this.S{nr}*x, t) * this.S{nr};
            
            % the multiplication by S{nr} from the right is due to the
            % inner derivative (chain rule).
            % this is identity if no projection occured.
            % If this function has been projected, the argument passed to
            % it is Vz, so the component set gradient is right-multiplied
            % by the components of V (which are in S in this case, see
            % setPointSet).
        end
        
        function J = evaluateJacobianSet(this, nr, x, t)
            % Returns the jacobian entries of the point set that have been
            % specified using setPointSet's argument jpd.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: A state space column vector `\vX` @type colvec<double>
            % t: The corresponding time `t` @type double
            %
            % Return values:
            % J: The values of the given points (linear indexing) in the
            % Jacobian in each row.
            %
            % See also: setPointSet
            
            J = this.evaluateComponentPartialDerivatives(this.PointSets{nr},...
                this.jend{nr}, this.jrow{nr}, this.deriv{nr}, ...
                this.jself{nr}, this.S{nr}*x, t, this.dfxsel{nr});
            % The deriv{nr} contains only derivative indices for non-zero
            % jacobian elements (determined in setPointSet using the
            % JSparsityPattern). Thus, the return values of
            % evaluateComponentPartialDerivatives only are of the size of
            % the actually non-zero jacobian values.
            % The T matrix is a sparse transformation, that places those
            % values into the correct spots with all other values that have
            % been demanded to be evaluated but are actually zero.
            %
            % See setPointSet
            J = this.T{nr} * J;
        end
        
        function J = evaluateJacobianSetMulti(this, nr, x, t, mu)
            % Returns the jacobian entries at multiple locations/times/parameters of the point set that have been
            % specified using setPointSet's argument jpd.
            %
            % Multi-argument version of evaluateJacobianSet.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
            % matrix<double>
            % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
            %
            % Return values:
            % J: The values of the given points (linear indexing) in the Jacobian in each row,
            % with as many columns as `\vX` had.
            %
            % See also: setPointSet
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
            J = zeros(length(this.deriv{nr}), m);
            for idx = 1:size(x,2)
                if ~singlemu
                    this.prepareSimulation(mu(:, idx));
                end
                J(:,idx) = this.evaluateComponentPartialDerivatives(this.PointSets{nr},...
                this.jend{nr}, this.jrow{nr}, this.deriv{nr}, ...
                this.jself{nr}, this.S{nr}*x(:,idx), t(idx), this.dfxsel{nr});
            end
            
            % The deriv{nr} contains only derivative indices for non-zero
            % jacobian elements (determined in setPointSet using the
            % JSparsityPattern). Thus, the return values of
            % evaluateComponentPartialDerivatives only are of the size of
            % the actually non-zero jacobian values.
            % The T matrix is a sparse transformation, that places those
            % values into the correct spots with all other values that have
            % been demanded to be evaluated but are actually zero.
            %
            % See setPointSet
            J = this.T{nr} * J;
        end  
               
        function setPointSet(this, nr, pts, jpd)
            % Parameters:
            % pts: A row vector with the desired points @type rowvec<integer>
            % jpd: ("Jacobian Partial Derivatives") A cell array of size
            % equal to the number of points. Each cell contains the indices
            % for which the partial derivatives of the corresponding
            % component function will be computed when calling
            % evaluateJacobianSet. @type cell
            
            if nr > 4
                warning('KerMor:Devel',['Point set numbering seems to ' ...
                    'be used in a different context. 1-4 are safe to work '...
                    'with at this stage, make sure nothing gets overridden.']);
            end
            % Ensure row vectors
            if size(pts,1) > 1
                pts = reshape(pts,1,[]);
            end
            if length(unique(pts)) ~= length(pts)
                error('Points have to be unique.');
            end
            
            this.pts{nr} = pts;
            
            jr = []; js = logical.empty(0,1);
            je = zeros(1,length(pts));
            SP = this.JSparsityPattern;
            if ~isempty(SP)
                % Jacobian extras
                if nargin == 4
                    deri = [];
                    full_mapping = [];
                    requested_len = 0;
                    dfx_sel = sparse(length(pts),0);
                end
                for i=1:length(pts)
                    sprow = SP(pts(i),:);
                    inew = find(sprow);
                    if isempty(inew)
                        inew = 1;
                        warning('setPointSet: Found all-zero Jacobian row for point %d. Using first state-space dof as "dependency"',pts(i));
                    end
                    jr = [jr inew];%#ok
                    js = [js inew == pts(i)];%#ok
                    je(i) = length(jr);
                    % Jacobian extras
                    if nargin == 4
                        des_der = jpd{i};
                        % Select elements of sparsity pattern that have been
                        % requested to be evaluated on evaluateJacobianSet
                        full_mapping = [full_mapping sprow(des_der)];%#ok
                        % deri contains the accumulated indices of 
                        % jacobian entries that are nonzero and thus make sense
                        % to evaluate.
                        [~, dpos, matchidx] = intersect(des_der, inew);
                        % Take the match indices, but sort them in the order
                        % the elements occur in jpd{i} to maintain correct
                        % ordering.
                        [~, sidx] = sort(dpos);
                        pos = reshape(matchidx(sidx),1,[]);

    %                     pos1 = jpd{i}(deriv_elem);
                        if ~isempty(pos)
                            if i==1
                                offs = 0;
                            else
                                offs = je(i-1);
                            end    
                            deri = [deri pos+offs];%#ok
                        end
                        % Sum up the total length for later transformation
                        requested_len = requested_len + length(des_der);

                        % Augment dfx selection matrix (only needed for default
                        % finite differences implementation)
                        dfx_sel(i,(end+1):(end+length(inew))) = 1;%#ok
                    end
                end
            else
                warning('Empty JSparsityPattern. DEIM approximation will not use any state space arguments.');
            end
            this.jrow{nr} = jr;
            this.jself{nr} = js;
            this.jend{nr} = je;
           
            % Compose x-entry selection matrix
            len = je(end);
            sel = jr(1:len);
            if ~isempty(this.V)
                this.S{nr} = this.V(sel,:);
            else
                this.S{nr} = sparse(1:len,sel,ones(len,1),len,this.xDim);
            end
            
            % Extras for jacobian evaluations
            this.deriv{nr} = [];
            this.T{nr} = [];
            if nargin == 4
                this.deriv{nr} = deri';
                % Find all the indices of the full mapping record that are
                % nonzero
                at = find(full_mapping);
                % Create a sparse transformation matrix, that assigns the
                % possibly nonzero jacobian entries returned by
                % evaluateComponentPartialDerivatives to their correct spot
                % according to the initially demanded derivatives (this
                % inserts zeros everywhere). If no zero jacobian entries
                % have been specified in setPointSet this is the identity
                % matrix.
                this.T{nr} = sparse(at,1:length(at),ones(length(at),1),...
                    requested_len,length(at));
                % Store convenience dfx selection matrix
                this.dfxsel{nr} = logical(dfx_sel);
            end
        end
        
        function target = project(this, V, W, target)
            if nargin < 4
                target = this.clone;
            end
            target = project@dscomponents.ACoreFun(this, V, W, target);
            % In case the point sets are not created newly: update the
            % current selection matrices
            for i=1:length(this.S)
                if ~isempty(this.S{i})
                    target.S{i} = this.S{i}*V;
                end
            end
        end
        
        function copy = clone(this, copy)
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.pts = this.pts;
            copy.S = this.S;
            copy.jrow = this.jrow;
            copy.jend = this.jend;
            copy.jself = this.jself;
            copy.deriv = this.deriv;
            copy.dfxsel = this.dfxsel;
            copy.T = this.T;
        end
        
        function psets = get.PointSets(this)
            psets = this.pts;
        end
        
        function res = test_ComponentEvalMatch(this, xsize)
            % Tests if the local implementation of evaluateComponents matches the full
            % evaluation
            
            if ~isempty(this.V)
                error('Cannot run test for projected ACompEvalCoreFuns.');
            end
            
            if nargin < 2
                xsize = 5;
            end
            x = rand(this.xDim,xsize);
            mu = rand(20,xsize); % simply assume param dim<=20
            t = rand(1,xsize);
            fx = this.evaluateMulti(x, t, mu);
            oldpts = [];
            if ~isempty(this.PointSets)
                oldpts = this.PointSets{1};
            end
            
            nsets = 3;
            s = RandStream('mt19937ar','Seed',2);
            % Limit set sizes to 5000
            %setsizes = s.randi(min(size(fx,1),5000), nsets, 1); @temp
            setsizes = s.randi(200, nsets, 1);
            sets = cell(nsets,1);
            for i = 1:length(setsizes)
                sets{i} = unique(s.randi(size(fx,1),1,setsizes(i)));
            end
            % Add an extra set with full size (only for functions with dim less than 10000)
            if size(fx,1) <= 10000
                sets{end+1} = 1:size(fx,1);
            end
            res = true;
	    
            for idx=1:length(sets)
                set = sets{idx};
                this.setPointSet(1, set);
                fxc = this.evaluateComponentSetMulti(1, x, t, mu);
                tmp = fx(set,:);
                err = abs((tmp-fxc)./tmp);
                err = err(tmp ~= 0);
                d = max(err);
                lres = isempty(d) || d < 1e-3;
                if ~lres || d > 10*sqrt(eps)
                    fprintf(2,'Warning! ACompEvalCoreFun evaluation test: Max rel. diff. for component wise vs. full evaluation: %e\n',full(d));
                end
                res = res && lres;
            end
            if ~isempty(oldpts)
                this.setPointSet(1,oldpts);
            end
        end
    end
    
    methods(Access=protected)
        function dfx = evaluateComponentGradientsAt(this, pts, ends, idx, self, x, t)
            % Default implementation of gradient computation via finite
            % differences.
            %
            % Override in subclasses for more performance if direct
            % derivative information is available.
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            %
            % Return values:
            % dfx: A length(pts) x size(x,1) matrix containing the
            % non-zero elements of the gradients of each component function
            % specified by pts. @type matrix<double>
            dt = sqrt(eps);
            d = size(x,1);
            X = repmat(x,1,d); T = repmat(t,1,d); MU = repmat(this.mu,1,d);
            I = speye(d,d)*dt;
            hlp = this.evaluateComponents(pts, ends, idx, self, x, t);
            dfx = (this.evaluateComponentsMulti(pts, ends, idx, self, X+I, T, MU) ...
                - repmat(hlp,1,d))/dt;
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, dfxsel)
            % Computes specified partial derivatives of `f` of the components given by pts and
            % the selected partial derivatives by dfxsel.
            %
            % Override in subclasses for optimized performance
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % deriv: The indices within `\vx` that derivatives are required for.
            % @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            % dfxsel: A derivative selection matrix. Contains the mapping for each row of x to
            % the output points pts. As deriv might contain less than 'size(x,1)' values, use
            % 'dfxsel(:,deriv)' to select the mapping for the actually computed derivatives.
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            %
            % See also: setPointSet
            dt = sqrt(eps);
            d = length(deriv);
            xd = size(x,2);
            
            X = repmat(x,1,d); T = repmat(t,1,d); %#ok<*PROP>
            I = sparse(deriv,1:d,ones(1,d),size(x,1),d)*dt;
            dfx = (this.evaluateComponentsMulti(pts, ends, idx, self, X+I, T, this.mu) ...
                - this.evaluateComponentsMulti(pts, ends, idx, self, X, T, this.mu))/dt;
            dfx = dfx(repmat(dfxsel(:,deriv),1,xd));
        end
        
        function dfx = evaluateComponentPartialDerivativesMulti(this, pts, ends, idx, deriv, self, x, t, mu, dfxsel)
            % Multi-argument evaluation method for partial derivatives.
            % Not used so far in KerMor, this is "legacy code" to keep
            % around if needed at any stage as default finite
            % difference-implementation.
            dt = sqrt(eps);
            d = length(deriv);
            xd = size(x,2);
            el = reshape(repmat(1:xd,d,1),1,[]);
            X = x(:,el); T = t(el); MU = mu(:,el);
            I = repmat(sparse(deriv,1:d,ones(1,d),size(x,1),d)*dt,1,xd);
            dfx = (this.evaluateComponentsMulti(pts, ends, idx, self, X+I, T, MU) ...
                - this.evaluateComponentsMulti(pts, ends, idx, self, X, T, MU))/dt;
            dfx = dfx(repmat(dfxsel(:,deriv),1,xd));
            dfx = reshape(dfx,[],xd);
        end

        function fx = evaluateComponentsMulti(this, pts, ends, idx, self, x, t, mu)
            % @todo improve performance!
            fx = zeros(length(pts),size(x,2));
            for k=1:size(x,2)
                this.prepareSimulation(mu(:,k));
                fx(:,k) = this.evaluateComponents(pts, ends, idx, self, x(:,k), t(k));
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % This is the template method that actually evaluates the components at given points
        % and values.
        %
        % @attention This method must be able to handle vector-arguments
        % for `\vx,t,\vmu`!
        %
        % Parameters:
        % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
        % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
        % indicates an input value relevant for the `i`-th point evaluation, i.e.
        % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
        % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
        % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
        % self: The positions in the `\vx` vector that correspond to the `i`-th output
        % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
        % @type rowvec<integer>
        % x: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
        % matrix<double>
        % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
        %
        % Return values:
        % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
        % many columns as `\vX` had.
        fx = evaluateComponents(this, pts, ends, idx, self, x, t, mu);
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if nargin == 2
                obj.pts = from.pts;
                obj.S = from.S;
                obj.jrow = from.jrow;
                obj.jend = from.jend;
                obj.jself = from.jself;
                obj.deriv = from.deriv;
                obj.dfxsel = from.dfxsel;
                obj.T = from.T;
                obj = loadobj@dscomponents.ACoreFun(obj, from);
            else
                obj = loadobj@dscomponents.ACoreFun(obj);
            end
        end
    end
end