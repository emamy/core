classdef BellFunction < kernels.ARBFKernel
    %BELLFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @todo export different estimator strategies into external classes =>
    % simulation speedup, separation of concerns..
    %
    % @todo investigate why the newton iteration sometimes exceeds max
    % iteration limit (example: intermittently test_LinearModelParams is
    % such a case)
    %
    % @docupdate Properties and class description
    %
    % @change{0,4,dw,2011-06-07} Moved the ModifiedNewton methods from
    % error.lipfun.ImprovedLocalSecantLipschitz to this class as they are more appropriate here.
    %
    % @change{0,4,dw,2011-05-27} Changed `x_0` to `r_0` and `x_R` to `r_m` as adopted in the WH10
    % Paper.
    %
    % @change{0,4,dw,2011-05-20} Removed the getLipschitzFunction method as it causes LARGE overhead
    % when being called very often. Instead, the estimator always uses the improved estimation
    % procedure.
    %
    % @change{0,4,dw,2011-05-19}
    % - Removed the PenaltyFactor property as no longer needed.
    % - Modified the newton iteration penalization such that the
    % 2nd degree polynomials extend the objective function both in value and derivative. This
    % ensures propert continuation of the objective function into the penalized area with respect to
    % the current kernel configuration (i.e. kernels.GaussKernel: small Gamma property causes high
    % derivatives and thus large estimations.
    %
    % @new{0,3,dw,2011-04-06} Added a new property
    % kernels.BellFunction.MaxNewtonIterations that avoids computations to
    % hang if the newton iteration does not come to a hold. An error will
    % be thrown as finding the correct minima is necessary.
    
    properties(SetObservable, Dependent)
        % Point of maximum first derivative on scalar evaluation.
        %
        % @propclass{critical} This value is essential for any bell function.
        r0;
    end
    
    properties(SetObservable)
        % Error tolerance for modified newton iteration
        %
        % @propclass{optional} The default precision is sufficient for all cases encountered so far.
        %
        % @default 1e-7
        NewtonTolerance = 1e-7;
        
        % Hard break for iteration count of modified newton algorithm
        %
        % @propclass{alglimit} Emergency break if iteration does not converge.
        %
        % @default 5000
        MaxNewtonIterations = 5000;
    end
    
    properties(SetAccess=private)
        % The maximum ("right") value for any `r_s`.
        rm;
    end
    
    properties(Access=private, Transient)
        p;
        sc;
        d1;
        d2;
    end
    
    properties(Access=private)
        fr0;
    end
    
    methods
        
        function this = BellFunction
            this = this@kernels.ARBFKernel;
            this.registerProps('r0','NewtonTolerance','MaxNewtonIterations');
        end
        
        function c = getGlobalLipschitz(this)
            % Computes the absolute value of the first derivative at x0
            % Implements the template method from BaseKernel.
            c = abs(this.evaluateD1(this.r0));
        end
        
        function copy = clone(this, copy)
            copy = clone@kernels.ARBFKernel(this, copy);
            copy.NewtonTolerance = this.NewtonTolerance;
            copy.MaxNewtonIterations = this.MaxNewtonIterations;
            copy.fr0 = this.fr0;
            copy.p = this.p;
            copy.sc = this.sc;
            copy.d1 = this.d1;
            copy.d2 = this.d2;
            copy.rm = this.rm;
            copy.r0 = this.r0;
        end
        
        function rtmp = ModifiedNewton(this, rstart, s)
            rtmp = rstart+2*this.NewtonTolerance;
            r = rstart;
            
            % Add eps at division as for zero nominator the expression is zero. (no error made)
            n = zeros(3,length(s));
            n(1,:) = this.d1(1) - (this.sc(1)-this.evaluateScalar(s))./(eps-s);
            n(2,:) = this.d1(2) - (this.sc(2)-this.evaluateScalar(s))./(this.p(2)-s+eps);
            n(3,:) = this.d1(3) - (this.sc(3)-this.evaluateScalar(s))./(this.p(3)-s+eps);
            
            % Lim x->y n'(x) = \phi''(x) / 2 (De L'hospital)
            dn = zeros(3,length(s));
            dn(1,:) = this.d2(1) - n(1,:)./-s;
            dn(1,isnan(dn(1,:))) = this.evaluateD2(0)/2;
            dn(2,:) = this.d2(2) - n(2,:)./(this.p(2)-s);
            dn(2,isinf(dn(2,:))) = 0; % == ddf(x0)/2;
            dn(3,:) = this.d2(3) - n(3,:)./(this.p(3)-s);
            dn(3,isinf(dn(3,:))) = this.d2(3)/2;
            
            cnt = 0;
            
            %conv = 1:length(x);
            %finished = zeros(size(x));
            %while ~isempty(conv) && cnt < this.MaxNewtonIterations
            while any(abs(rtmp-r) > this.NewtonTolerance) && cnt < this.MaxNewtonIterations
                
                % Experimental setting: thinning out the finished values does not improve the
                % overall performance.
                % Find indices of finished values
                %                 fidx = abs(xtmp-x) <= this.NewtonTolerance;
                %                 if any(fidx)
                %                     % Store finished values
                %                     finished(conv(fidx)) = xtmp(fidx);
                %                     conv(fidx) = [];
                %                     x(fidx) = []; xtmp(fidx) = []; y(fidx) = [];
                %                     n0(fidx) = []; dn0(fidx) = [];
                %                     nx0(fidx) = []; dnx0(fidx) = [];
                %                     nxr(fidx) = []; dnxr(fidx) = [];
                %                 end
                
                [g,dg] = this.optFun(r,s,n,dn);
                
                rtmp = r;
                r = r - g./dg;
                cnt = cnt + 1;
                
                %if cnt > this.MaxNewtonIterations*.5
                %showNewton;
                %end
            end
            if cnt == this.MaxNewtonIterations
                error('Bellfunction->ModifiedNewton: Max iterations of %d reached',this.MaxNewtonIterations);
            end
            
            %             function showNewton %#ok<*DEFNU>
            %                 % Fix an identical y for all iterations
            %                 si = 200;
            %                 for i=1:5
            %
            %                     X = linspace(.5*min(min([x; xtmp],[],1)),2*max(max([x; xtmp],[],1)),si);
            %                     gs = zeros(length(y),si);
            %                     for idx = 1:length(y)
            %                         gs(idx,:) = error.lipfun.ImprovedLocalSecantLipschitz.optFun(X,repmat(y(idx),size(X)),f,df,ddf,this.x0);
            %                     end
            %
            %                     [g,dg] = error.lipfun.ImprovedLocalSecantLipschitz.optFun(x,y,f,df,ddf,this.x0);
            %                     xtmp = x;
            %                     x = x - g./dg;
            %                     %plot(xs,gs,'r');%,'LineWidth',2
            %                     figure(1);
            %                     subplot(1,2,1);
            %                     plot(X,gs,'b',X,0,'black-');
            %                     hold on;
            %                     plot([xtmp; x],[g; zeros(size(g))],x,0,'rx',xtmp,g,'rx');
            %                     hold off;
            %                     axis tight;
            %
            %                     subplot(1,2,2);
            %                     X = linspace(0-this.xR,2*this.xR,si);
            %                     gs = zeros(length(y),si);
            %                     for idx = 1:length(y)
            %                         gs(idx,:) = kernels.BellFunction.optFun(X,repmat(y(idx),size(X)),f,df,ddf,this.x0,this.PenaltyFactor);
            %                     end
            %                     plot(X,gs);
            %                     hold on;
            %                     plot(X,0,'black-');
            %                     hold off;
            %                     pause;
            %                 end
            %             end
        end
        
        function set.r0(this, value)
            if ~isreal(value) || value < 0 || ~isscalar(value)
                error('r0 must be a scalar greater than zero.');
            end
            this.fr0 = value;
            this.setConstants;
        end
        
        function value = get.r0(this)
            value = this.fr0;
        end
    end
    
    methods(Access=private)
        function setConstants(this)
            % Sets constants depending only on the current `r_0` value.
            
            % Set outer bound for 'r_s' positions
            this.rm = this.evaluateScalar(0)*this.r0 /...
                (this.evaluateScalar(0)-this.evaluateScalar(this.r0));
            % Set values needed for Newton iteration
            q = [0 this.r0 this.rm];
            this.sc = this.evaluateScalar(q);
            this.d1 = this.evaluateD1(q);
            this.d2 = this.evaluateD2(q);
            this.p = q;
        end
        
    end
    
    methods(Access=?demos.BellFunctions)
        
        function [g,dg,pi,pl,pr] = optFun(this,r,s,n,dn)
            g = zeros(size(r));
            d = g; c = d; dg = g;
            
            a = s < this.r0;
            pi = a & r < this.r0 | ~a & r > this.r0;
            pr = a & r > this.rm;
            pl = ~a & r < 0;
            std = ~(pi | pl | pr);
            
            % Standard case
            rs = r(std); ss = s(std);
            g(std) = this.evaluateD1(rs) - (this.evaluateScalar(rs)-this.evaluateScalar(ss))./(rs-ss);
            g(isnan(g)) = 0;
            dg(std) = this.evaluateD2(rs) - g(std)./(rs-ss);
            ninf = isinf(dg);
            dg(ninf) = this.evaluateD2(rs(ninf))/2;
            
            % Penalty case
            penalty = ~std;
            c(pi) = dn(2,pi).^2 ./ (4*n(2,pi));
            d(pi) = 2*n(2,pi)./dn(2,pi) - this.r0;
            c(pl) = dn(1,pl).^2 ./ (4*n(1,pl));
            d(pl) = 2*n(1,pl)./dn(1,pl);
            c(pr) = dn(3,pr).^2 ./ (4*n(3,pr));
            d(pr) = 2*n(3,pr)./dn(3,pr) - this.rm;
            
            g(penalty) = c(penalty) .* (r(penalty) + d(penalty)).^2;
            dg(penalty) = 2 * c(penalty) .* (r(penalty) + d(penalty));
            % Avoid exact matches!
            g(dg == 0) = 0;
            dg(dg == 0) = 1;
        end
    end
    
    methods(Abstract)
        % Method for first derivative evaluation
        dr = evaluateD1(r);
        
        % Method for second derivative evaluation
        ddr = evaluateD2(r);
    end
    
    methods(Static)
        
        
        %         function showNewton(x, x0, ybase, f, df, ddf)
        %             % Fix an identical y for all iterations
        %             y = repmat(ybase,1,size(x,2));
        %             gs = kernels.BellFunction.optFun(x,y,f,df,ddf,x0);
        %             xs = x;
        %             for i=1:5
        %                 [g,dg] = error.lipfun.ImprovedLocalSecantLipschitz.optFun(x,y,f,df,ddf,x0);
        %                 xtmp = x;
        %                 x = x - g./dg;
        %                 plot(xs,gs,'r');%,'LineWidth',2
        %                 hold on;
        %                 plot([xtmp; x],[g; zeros(size(g))],x,0,'rx');
        %                 hold off;
        %                 title(sprintf('ybase=%f',ybase));
        %                 %axis([min(xs) max(xs) min(gs)*1.5 max(gs)*1.5]);
        %                 pause;
        %             end
        %         end
        
    end
    
    methods(Static,Access=protected)
        function this = loadobj(this, from)
            % As the constant properties are transient, they have to be re-computed upon loading.
            %
            % Any subclasser MUST call this superclasses loadobj function explicitly!
            if nargin > 1
                % NOTE: r0 must be set via setter in subclass at some stage.
                if isfield(from,'fr0') && ~isempty(from.fr0)
                    this.fr0 = from.fr0;
                end
                this.MaxNewtonIterations = from.MaxNewtonIterations;
                this.NewtonTolerance = from.NewtonTolerance;
                this = loadobj@KerMorObject(this, from);
            elseif ~isa(this, 'kernels.BellFunction')
                error('Object passed is not a kernels.BellFunction instance and no init struct is given.');
            end
            if isempty(this.fr0)
                error('Internal r0 field not persisted yet and not initialized in subclass loadobj method.');
            end
            this.setConstants;
        end
    end
    
end

