classdef BellFunction < kernels.BaseKernel & kernels.IRotationInvariant
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
    % @new{0,3,dw,2011-04-06} Added a new property
    % kernels.BellFunction.MaxNewtonIterations that avoids computations to
    % hang if the newton iteration does not come to a hold. An error will
    % be thrown as finding the correct minima is necessary.
    
    properties(SetObservable)
        % Point of maximum first derivative on scalar evaluation.
        %
        % @propclass{critical} This value is essential for any bell function.
        x0;
        
        % Penalty factor for modified newton iteration.
        %
        % @propclass{optional} Value is set to 2 by developer (empirically).
        %
        % @default 2
        PenaltyFactor = 2;
        
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
    
    properties(SetAccess=private, Dependent)
        % The maximum ("right") value for any `x_y`.
        xR;
    end
    
    properties(Access=private)
        oldxfeat = [];
    end
    
    properties(Access=private, Transient)
        priv_xr = [];
    end
    
    methods
        
        function this = BellFunction
            this = this@kernels.BaseKernel;
            this.registerProps('x0','PenaltyFactor','NewtonTolerance','MaxNewtonIterations');
        end
        
        function fcn = getLipschitzFunction(this)
            % Overrides the method from the base class kernels.BaseKernel
            % as for this type of kernels a better local Lipschitz constant
            % estimation is available.
            %
            % See also: kernels.BaseKernel#getLipschitzFunction
            fcn = @this.getImprovedLocalSecantLipschitz;
        end
        
        function c = getGlobalLipschitz(this)
            % Computes the absolute value of the first derivative at x0
            % Implements the template method from BaseKernel.
            c = abs(this.evaluateD1(this.x0));
        end
        
        function ci = getLocalGradientLipschitz(this, di, C, t, mu)
            
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                ci = ones(size(di))*abs(this.evaluateD1(this.x0));
            else
                case1 = di - C - this.x0 > 0;
                case2 = di + C - this.x0 < 0;
                ci(case1) = abs(this.evaluateD1(di(case1)-C));
                ci(case2) = abs(this.evaluateD1(di(case2)+C));
                ci(~case1 & ~case2) = abs(this.evaluateD1(this.x0));
            end
        end
        
        function ci = getLocalSecantLipschitz(this, di, C, t, mu)
            
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                ci = ones(size(di))*abs(this.evaluateD1(this.x0));
            else
                case1 = di - C - this.x0 > 0;
                case2 = di + C - this.x0 < 0;
                
                % If C is too small we just take the derivative at this
                % point.
                if C < sqrt(eps)
                    ci(case1 | case2) = abs(this.evaluateD1(di(case1 | case2)));
                else
                    ci(case1) = (this.evaluateScalar(di(case1)-C) - this.evaluateScalar(di(case1))) / C;
                    ci(case2) = (this.evaluateScalar(di(case2)) - this.evaluateScalar(di(case2)+C)) / C;
                end
                ci(~case1 & ~case2) = abs(this.evaluateD1(this.x0));
            end
        end
        
        function ci = getImprovedLocalSecantLipschitz(this, di, C, t, mu)
            if isempty(this.oldxfeat) || any(isnan(this.oldxfeat)) || size(this.oldxfeat,2) ~= size(di,2)
                this.oldxfeat = this.x0+.2*sign(this.x0-di);
            end
            
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                xfeat = this.ModifiedNewton(this.oldxfeat,di);
                ci = abs(this.evaluateD1(xfeat));
                this.oldxfeat = xfeat;
            else
                xfeat = ones(size(di))*this.x0;
                update = abs(di-this.x0) < C;
                % Choose suitable starting conditions if no old xfeat vectors
                % are available
                if any(update)
                    xfeat(update) = this.ModifiedNewton(this.oldxfeat(update),di(update));
                    this.oldxfeat(update) = xfeat(update);
                end
                left = di + C - xfeat < 0;
                right = di - C - xfeat > 0;
                center = ~left & ~right;
                % If C is too small we just take the derivative at this
                % spot
                if C < sqrt(eps)
                    ci(left | right) = abs(this.evaluateD1(di(left | right)));
                else
                    ci(left) = (this.evaluateScalar((di(left))) - this.evaluateScalar((di(left)+C))) / C;
                    ci(right) = (this.evaluateScalar((di(right)-C)) - this.evaluateScalar((di(right)))) / C;
                end
                ci(center) = abs(this.evaluateD1(xfeat(center)));
            end
        end
        
        function set.x0(this, value)
            if ~isposrealscalar(value) 
                error('x0 must be a scalar greater than zero.');
            end
            this.x0 = value;
            % Ensure to reset the current newton iteration starting point
            % to avoid having endless loops in newton iterations.
            this.oldxfeat = [];%#ok
        end
        
        function set.PenaltyFactor(this, value)
            this.PenaltyFactor = value;
            keyboard;
        end
        
        function value = get.xR(this)
            if isempty(this.priv_xr)
                this.priv_xr = this.evaluateScalar(0)*this.x0 /...
                    (this.evaluateScalar(0)-this.evaluateScalar(this.x0));
            end
            value = this.priv_xr;
        end
    end
    
    methods(Access=private)
        
        function xtmp = ModifiedNewton(this, xstart, y)
            xtmp = xstart+2*this.NewtonTolerance;
            x = xstart;
            f = @(x)this.evaluateScalar(x);
            df = @(x)this.evaluateD1(x);
            ddf = @(x)this.evaluateD2(x);
            
            n0 = df(0) - (f(0)-f(y))./-y;
            nx0 = df(this.x0) - (f(this.x0)-f(y))./(this.x0-y);
            % Numerical fix: Sometimes y equals x0, in that case nx0 = 0.
            nx0(isnan(nx0)) = 0;
            nxr = df(this.xR) - (f(this.xR)-f(y))./(this.xR-y);
            
            p = this.PenaltyFactor;
            cnt = 0;
            while any(abs(xtmp-x) > this.NewtonTolerance) && cnt < this.MaxNewtonIterations
                
                [g,dg] = optFun(x,y);
                
                xtmp = x;
                x = x - g./dg;
                cnt = cnt + 1;
                
                %showNewton;
            end
            if cnt == this.MaxNewtonIterations
                error('Bellfunction->ModifiedNewton: Max iterations of %d reached',this.MaxNewtonIterations);
            end
            
%             function showNewton
%                 % Fix an identical y for all iterations
%                 si = 200;
%                 for i=1:5
%                     
%                     X = linspace(.5*min(min([x; xtmp],[],1)),2*max(max([x; xtmp],[],1)),si);
%                     gs = zeros(length(y),si);
%                     for idx = 1:length(y)
%                         gs(idx,:) = kernels.BellFunction.optFun(X,repmat(y(idx),size(X)),f,df,ddf,this.x0,this.PenaltyFactor);
%                     end
%                     
%                     [g,dg] = kernels.BellFunction.optFun(x,y,f,df,ddf,this.x0,this.PenaltyFactor);
%                     xtmp = x;
%                     x = x - g./dg;
%                     %plot(xs,gs,'r');%,'LineWidth',2
%                     plot(X,gs,'b',X,0,'black-');
%                     hold on;
%                     plot([xtmp; x],[g; zeros(size(g))],x,0,'rx',xtmp,g,'rx');
%                     hold off;
%                     axis tight;
%                     pause;
%                 end
%             end
            
            function [g,dg] = optFun(x,y)
                g = zeros(size(x));
                dg = g;
                
                a = y <= this.x0;
                std = (a & x < this.xR & x > this.x0) | (~a & 0 < x & x < this.x0);
                p1 = a & x <= this.x0 | ~a & x > this.x0;
                p4 = a & x >= this.xR;
                p2 = ~a & x <= 0;
                
                % Standard case
                xs = x(std); ys = y(std);
                g(std) = df(xs) - (f(xs)-f(ys))./(xs-ys);
                g(isnan(g)) = 0;
                dg(std) = ddf(xs) - g(std)./(xs-ys);
                if any(isnan(dg))
                    dg(isnan(dg)) = eps%#ok
                end
                
                if any(p1)
                    g(p1) = nx0(p1) - p*(x(p1)-this.x0).^2;
                    dg(p1) = -2*p*(x(p1)-this.x0);
                end
                if any(p2)
                    g(p2) = n0(p2) + p*x(p2).^2;
                    dg(p2) = 2*p*x(p2);
                end
                if any(p4)
                    g(p4) = nxr(p4) + p*(x(p4)-this.xR).^2;
                    dg(p4) = 2*p*(x(p4)-this.xR);
                end
            end
            
        end
        
    end
    
    methods(Abstract)
        % Method for first derivative evaluation
        dx = evaluateD1(x);
        
        % Method for second derivative evaluation
        ddx = evaluateD2(x);
    end
    
    methods(Static)
        function xFeatDemo(Gamma)
            if nargin == 0
                Gamma = 2;
            end
            
            bell = kernels.GaussKernel(Gamma);
            bell.NewtonTolerance = 1e-3;
            h = figure(1);
            
            maxX = Gamma*3;
            x = linspace(0,maxX,70);
            
            %xstart = x0*(1+.5*sign(x0-x));
            xstart = repmat(bell.x0/2,1,size(x,2));
            
            xfeats = bell.ModifiedNewton(xstart,x);
            
            f = @bell.evaluateScalar;
            df = @bell.evaluateD1;
            ddf = @bell.evaluateD2;
            
            %kernels.BellFunction.showNewton(x,bell.x0,y,bell.PenaltyFactor,f,df,ddf);
            
            for k=1:length(x)
                y = x(k);
                % Get precomputed values
                dxf = df(xfeats);
                
                kernels.BellFunction.showNewton(x,bell.x0,y,bell.PenaltyFactor,f,df,ddf);
                
                % Get optimization function
                opt = kernels.BellFunction.optFun(x,repmat(y,1,size(x,2)),f,df,ddf,bell.x0,bell.PenaltyFactor);
                
                % Bound plot
                hlp = (maxX-y)*dxf(k); % also used in plot of tangent!
                opt(opt < hlp) = hlp;
                hlp2 = max(1.5,max(df(x)));
                opt(opt > hlp2) = hlp2;
                
                % Plot!
                plot(x,f(x),'r',x,opt,'b',x,abs(df(x)),'m',x,abs(dxf),'g',[0 y maxX],f(y) + [-y*dxf(k) 0 hlp],...
                    'g--',y,f(y),'r.',[xfeats(k) xfeats(k) y],[f(xfeats(k)) abs(dxf(k)) abs(dxf(k))],'blackx',x,0,'black');
                legend('Gaussian','Opt. target','Absolute Gaussian derivative','Max local derivatives','Max local gradient around x0','x0','xfeat (correlation)');
                pause;
            end
            
            close(h);
            
        end
        
        function [g,dg] = optFun(x,y,f,df,ddf,x0,p)
            g = zeros(size(x));
            dg = g;
            
            xr = f(0)*x0 / (f(0)-f(x0));
            n0 = df(0) - (f(0)-f(y))./-y;
            nx0 = df(x0) - (f(x0)-f(y))./(x0-y);
            nxr = df(xr) - (f(xr)-f(y))./(xr-y);
            
            a = y <= x0;
            std = (a & x < xr & x > x0) | (~a & 0 < x & x < x0);
            p1 = a & x <= x0 | ~a & x > x0;
            p4 = a & x >= xr;
            p2 = ~a & x <= 0;
            
            % Standard case
            xs = x(std); ys = y(std);
            g(std) = df(xs) - (f(xs)-f(ys))./(xs-ys);
            g(isnan(g)) = 0;
            dg(std) = ddf(xs) - g(std)./(xs-ys);
            if any(isnan(dg))
                dg(isnan(dg)) = eps%#ok
            end
            
            if any(p1)
                g(p1) = nx0(p1) - p*(x(p1)-x0).^2;
                dg(p1) = -2*p*(x(p1)-x0);
            end
            if any(p2)
                g(p2) = n0(p2) + p*x(p2).^2;
                dg(p2) = 2*p*x(p2);
            end
            if any(p4)
                g(p4) = nxr(p4) + p*(x(p4)-xr).^2;
                dg(p4) = 2*p*(x(p4)-xr);
            end
        end
        
        function showNewton(x, x0, ybase, p, f, df, ddf)
            % Fix an identical y for all iterations
            y = repmat(ybase,1,size(x,2));
            gs = kernels.BellFunction.optFun(x,y,f,df,ddf,x0,p);
            xs = x;
            for i=1:5
                [g,dg] = kernels.BellFunction.optFun(x,y,f,df,ddf,x0,p);
                xtmp = x;
                x = x - g./dg;
                plot(xs,gs,'r');%,'LineWidth',2
                hold on;
                plot([xtmp; x],[g; zeros(size(g))],x,0,'rx');
                hold off;
                title(sprintf('ybase=%f',ybase));
                %axis([min(xs) max(xs) min(gs)*1.5 max(gs)*1.5]);
                pause;
            end
        end
    end
end

