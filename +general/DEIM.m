classdef DEIM < KerMorObject & general.AProjectable & IReductionSummaryPlotProvider
    % DEIM: Implements the DEIM-Algorithm from \cite CS10
    %
    % Additionally allows to estimate the DEIM error a-posteriori by using the
    % next `M'` orders of the current DEIM approximation.
    %
    % @author Daniel Wirtz @date 2012-03-26
    %
    % @change{0,6,dw,2012-06-08}
    % - DEIM is now a KerMorObject
    % - Fixed triggering of updateOrderData after computeDEIM calls
    % - Completed clone method
    % - New event "OrderUpdated" to notify any components of changes to the
    % current DEIM orders
    %
    % @new{0,6,dw,2012-05-30} Added a custom implementation of the
    % getStateJacobian method
    %
    % @change{0,6,dw,2012-05-29} Besides many work-in-progress changes, now the
    % evaluate function is implemented directly. further, some matrices needed
    % for the DEIMEstimator have been included and are being computed.
    %
    % @new{0,6,dw,2012-03-26} Added this class.
    %
    % @todo think of exporting the jrow, jend, S properties to ACompEvalCoreFun
    % and precompute stuff there; projection happens at that stage
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetObservable)
        % The maximum order up to which the DEIM approximation should be computed.
        %
        % This corresponds to the maximum number `m` that can be choosen as approximation order.
        %
        % @propclass{important} The larger the maximal order, the better
        % approximation quality can be achieved, for paying the cost of
        % higher evaluation time.
        %
        % @type integer @default 40
        MaxOrder = 40;
    end
    
    properties(Dependent)
        % The actual order `M` for the current DEIM approximation.
        %
        % As second element, the order `M'` of the error estimation DEIM
        % approximation can be specified. If not given, `M'=0` is
        % automatically set.
        %
        % See also: MaxOrder
        %
        % @type rowvec<integer> @default []
        Order;
    end
    
    properties(SetAccess=private)
        % The singular values returned by the SVD decomposition to compute
        % the DEIM POD basis.
        %
        % This value is set when general.DEIM.computeDEIM is called.
        %
        % @type rowvec<double> @default []
        %
        % See also: computeDEIM
        SingularValues = [];
    end
    
    properties(GetAccess=protected, SetAccess=private)
        % The full approximation base
        %
        % @type matrix<double>
        u;
        
        % Interpolation points
        %
        % @type rowvec<double>
        %
        % See also: getInterpolationPoints
        pts;
    end
    
    properties(Access=private)
        % Vector with 2 entries for order of DEIM approximation and order
        % for error esimation
        %
        % @type rowvec<integer> @default []
        fOrder = [];
    end
    
    properties(SetAccess=protected)
        % The U matrix for the current Order.
        %
        % @type matrix<double>
        U;
        
        % If projection is applied, this contains the non-projected full
        %`d \times m` matrix `U_m(P_m^tU_m)^{-1}` for use in subclasses.
        %
        % @type matrix<double>
        U_nonproj;
        
        % Some matrices for M+M' error estimation
        %
        % @type matrix<double>
        Uerr1;
        Uerr2;
        
        M1;
        M2;
    end
    
    properties(SetAccess=private)
        % The function which DEIM is applied to
        %
        % Is a subclass of dscomponents.ACompEvalCoreFun
        %
        % @type dscomponents.ACompEvalCoreFun
        f;
        
        % The maximum residuals obtained along the magic points computation
        %
        % @type rowvec<double>
        Residuals;
    end
    
    methods
        function computeDEIM(this, f, fxi)
            % Implementation of the DEIM algorithm.
            %
            % Parameters:
            % f: the function DEIM is applied to @type dscomponents.ACompEvalCoreFun
            % fxi: the snapshot data @type matrix<double> or general.ABlockSVD
            
            this.f = f;
            if ~isa(this.f,'dscomponents.ACompEvalCoreFun');
                error('Cannot use DEIM with non ACompEvalCoreFun-implementing functions.');
            elseif ~f.test_ComponentEvalMatch
                error('Component evaluation does not match direct evaluation.');
            elseif ~f.test_ComponentEvalMatch(10)
                error('Component evaluation does not match direct evaluation for multi-arguments.');
            end
            
            if size(fxi,1) < this.MaxOrder
                fprintf(2,'Warning: MaxOrder %d of DEIM approximation is larger than actual f dimension %d. Setting MaxOrder=%d\n',...
                    this.MaxOrder,size(fxi,1),size(fxi,1));
                this.MaxOrder = size(fxi,1);
            end
            
            %% Generate u_1 ... u_m base
            p = general.POD;
            p.Mode = 'abs';
            p.Value = this.MaxOrder;
            
            %p.UseSVDS = size(fxi,1) > 10000;
            [Utmp, this.SingularValues] = p.computePOD(fxi);
            if isa(Utmp,'data.FileMatrix')
                Utmp = Utmp.toMemoryMatrix;
            end
            this.u = Utmp;
            
            if size(this.u,2) < this.MaxOrder
                fprintf('POD returned less (=%d) than MaxOrder (=%d) basis vectors. Setting MaxOrder=%d.\n',...
                    size(this.u,2),this.MaxOrder,size(this.u,2));
                this.MaxOrder = size(this.u,2);
            end
            
            % Get interpolation points
            this.pts = this.getInterpolationPoints(this.u);
            
            % Trigger computation of U matrix.
            % If no order has been set yet ("first computation"), use
            % predefined order. Otherwise, just ensure the order data
            % corresponds to the current u/pts.
            if isempty(this.fOrder)
                this.Order = ceil(this.MaxOrder/10);
            else
                this.updateOrderData;
            end
        end
        
        function fx = evaluate(this, x, t)
            fx = this.U * this.f.evaluateComponentSet(1, x, t);
        end
        
        function fx = evaluateMulti(this, x, t, mu)
            fx = this.U * this.f.evaluateComponentSetMulti(1, x, t, mu);
        end

        function J = getStateJacobian(this, x, t)
            hlp = this.U * this.f.evaluateComponentSetGradientsAt(1, x, t);
            if isempty(this.V)
                p = logical(this.f.JSparsityPattern);
                J = sparse(this.f.fDim,this.f.xDim);
                J(p) = hlp(p);
            else
                J = hlp;
            end
        end
        
        function err = getEstimatedError(this, x, t, mu)
            if this.fOrder(2) == 0
                error('No error estimation possible with zero ErrorOrder property');
            end
            err = this.Uerr1 * this.f.evaluateComponentSetMulti(1, x, t, mu) ...
                - this.Uerr2 * this.f.evaluateComponentSetMulti(2, x, t, mu);
        end
        
        function target = project(this, V, W, target)
            % Pojects instance according to the projection biorthogonal matrices `V,W`.
            %
            % See also: general.AProjectable
            %
            % Parameters:
            % V: The `V` matrix of the biorthogonal pair `V,W` @type
            % matrix<double>
            % W: The `W` matrix of the biorthogonal pair `V,W` @type
            % matrix<double>
            % target: Specify the subclasses projection target instance if
            % project is overridden in a subclass and the subclass has been
            % cloned. @type DEIM @default this
            
            if nargin < 4
                target = this.clone;
            end
            target = project@general.AProjectable(this, V, W, target);
            % Important: Project the component evaluable function, too!
            target.f = this.f.project(V, W);
            target.updateOrderData;
        end
        
        function copy = clone(this, copy)
            if nargin < 2
                copy = general.DEIM;
            end
            copy = clone@general.AProjectable(this, copy);
            % Clone associated f as different orders for the cloned object
            % lead to different PointSets for the component evaluable
            % function.
            copy.f = this.f.clone;
            copy.fOrder = this.fOrder;
            copy.u = this.u;
            copy.pts = this.pts;
            copy.U = this.U;
            copy.U_nonproj = this.U_nonproj;
            copy.MaxOrder = this.MaxOrder;
            copy.Uerr1 = this.Uerr1;
            copy.Uerr2 = this.Uerr2;
            copy.M1 = this.M1;
            copy.M2 = this.M2;
            copy.SingularValues = this.SingularValues;
        end
        
        function plotSummary(this, pm, context)
            if ~isempty(this.SingularValues)
                str = sprintf('%s: DEIM POD singular value decay, MaxOrder: %d',context,this.MaxOrder);
                h = pm.nextPlot('deim_singvals',...
                    str,'POD size','singular values');
                semilogy(h,this.SingularValues,'LineWidth',2);
            end
            if ~isempty(this.Residuals)
                str = sprintf('%s: maximum basis residual, MaxOrder: %d',context,this.MaxOrder);
                h = pm.nextPlot('max_residuals',...
                    str,'Point number','residual');
                semilogy(h,this.Residuals,'LineWidth',2);
            end
        end
    end
    
    methods(Access=private)
        function pts = getInterpolationPoints(this, u)
            % Computes the interpolation indices according to the DEIM
            % algorithm
            %
            % Parameters:
            % u: matrix, columns are the orthonormal basis vectors computed via POD from the snapshot data
            % @type matrix<double>
            
            n = size(u,1);
            m = size(u,2);
            pts = zeros(1, m);
            v = pts;
            [v(1), pts(1)] = max(abs(u(:,1)));
            P = zeros(n,1);
            P(pts(1)) = 1;
            for i=2:m
                c = (P'*u(:,1:(i-1))) \ (P'*u(:,i));
                [v(i), pts(i)] = max(abs(u(:,i) - u(:,1:(i-1))*c));
                P = sparse(pts(1:i),1:i,ones(i,1),n,i);
            end
            this.Residuals = v;
            if KerMor.App.Verbose > 2
                fprintf('DEIM interpolation points [%s] with values [%s]\n',...
                    Utils.implode(pts,' ','%d'),Utils.implode(v,' ','%2.2e'));
            end
        end
    end
    
    methods(Access=protected)
        function updateOrderData(this)
            % Update approximation order as specified in fOrder. As a
            % consequence some matrices have to be recalculated.
            
            if KerMor.App.Verbose > 3
                fprintf('general.DEIM.updateOrderData: Updating order data of DEIM (%s, #%s) to [%d %d]\n',class(this),this.ID,this.fOrder);
            end
            
            o = this.fOrder(1);
            n = size(this.u,1);
            
            P = sparse(this.pts(1:o),1:o,ones(o,1),n,o);
            this.U_nonproj = this.u(:,1:o) / (P'*this.u(:,1:o));
            
            % Set primary point set in ACompEvalCoreFun
            this.f.setPointSet(1, this.pts(1:o));
            
            om = this.fOrder(2);
            if om > 0
                % Use second point set in ACompEvalCoreFun
                this.f.setPointSet(2, this.pts(o+1:o+om));
                
                Perr = sparse(this.pts(o+1:o+om),1:om,ones(om,1),n,om);
                
                Um = this.u(:,1:o);
                Umd = this.u(:,o+1:o+om);
                A = P'*Um;
                B = P'*Umd;
                C = Perr'*Um;
                D = Perr'*Umd;
                E = C / A;
                F = D - E*B;
                this.Uerr2 = ((Um/A)*B - Umd) / F;
                this.Uerr1 = this.Uerr2 * E;
            else
                this.Uerr1 = [];
                this.Uerr2 = [];
            end
            
            if ~isempty(this.W)
                this.U = this.W'*this.U_nonproj;
                if om > 0
                    % Compute values for error estimator
                    this.M1 = this.Uerr1 + (Um - this.V*(this.W'*Um))/A;
                    this.M2 = this.Uerr2;
                    
                    % Project DEIM error estimation values
                    this.Uerr1 = this.W'*this.Uerr1;
                    this.Uerr2 = this.W'*this.Uerr2;
                end
            else
                this.U = this.U_nonproj;
            end
            
            % Fire order updated event
            notify(this, 'OrderUpdated', event.EventData);
        end
    end
    
    %% Getter & Setter
    methods
        function o = get.Order(this)
            o = this.fOrder;
        end
        
        function set.Order(this, value)
            if isempty(this.u) || isempty(this.pts)
                error('Cannot set DEIM order as approximateSystemFunction has not been called yet.');
            end
            if numel(value) ~= 2
                if numel(value) ~= 1
                    error('Order must be either a scalar or two element integer vector');
                else
                    value(2) = 0;
                end
            end
            if any(value) > size(this.u,2) || any(value) < 1
                error('Invalid Order/ErrOrder value. Allowed are integers in [1, %d]',size(this.u,2));
            elseif sum(value) > this.MaxOrder
                warning('KerMor:DEIM',sprintf(...
                    'Order (%d) and ErrOrder (%d) values may not exceed MaxOrder (%d). Using [%d 0]',...
                    value,this.MaxOrder,this.MaxOrder));%#ok
                value = [this.MaxOrder 0];
            end
            % Only re-compute error order matrices if an actual change
            % happenend
            if ~isequal(this.fOrder, value)
                this.fOrder = value;
                this.updateOrderData;
            end
        end
        
        function set.MaxOrder(this, value)
            if ~isposintscalar(value)
                error('MaxOrder has to be a positive integer scalar.');
            end
            this.MaxOrder = value;
        end
    end
    
    events
        % Gets fired whenever this DEIM instance has updated it's order
        % matrices.
        OrderUpdated;
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, varargin)
            obj = loadobj@general.AProjectable(obj, varargin{:});
            obj = loadobj@KerMorObject(obj, varargin{:});
        end
    end
end