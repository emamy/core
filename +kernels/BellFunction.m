classdef BellFunction < kernels.BaseKernel & kernels.IRotationInvariant
    %BELLFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0;
        
        PenaltyFactor = 2;
        
        NewtonTolerance = 1e-7;
    end
    
    properties(Access=private)
        oldxfeat = [];
    end
    
    %     properties(Dependent)
    %         xR;
    %     end
    
    methods
        
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
            if isempty(this.oldxfeat) || any(isnan(this.oldxfeat))
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
    end
    
    methods(Access=private)
        
        function xtmp = ModifiedNewton(this, xstart, y)
            xtmp = xstart+2*this.NewtonTolerance;
            x = xstart;
            f = @(x)this.evaluateScalar(x);
            df = @(x)this.evaluateD1(x);
            ddf = @(x)this.evaluateD2(x);
            
            xr = f(0)*this.x0 / (f(0)-f(this.x0));
            n0 = df(0) - (f(0)-f(y))./-y;
            nx0 = df(this.x0) - (f(this.x0)-f(y))./(this.x0-y);
            % Numerical fix: Sometimes y equals x0, in that case nx0 = 0.
            nx0(isnan(nx0)) = 0;
            nxr = df(xr) - (f(xr)-f(y))./(xr-y);
            
            p = this.PenaltyFactor;
            while any(abs(xtmp-x) > this.NewtonTolerance)
                
                [g,dg] = optFun(x,y);
                
                xtmp = x;
                x = x - g./dg;
            end
            
            function [g,dg] = optFun(x,y)
                g = zeros(size(x));
                dg = g;
                
                a = y <= this.x0;
                std = (a & x < xr & x > this.x0) | (~a & 0 < x & x < this.x0);
                p1 = a & x <= this.x0 | ~a & x > this.x0;
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
                    g(p1) = nx0(p1) - p*(x(p1)-this.x0).^2;
                    dg(p1) = -2*p*(x(p1)-this.x0);
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
            
            %showNewton(x,x0,c,f,df,ddf);
            
            for k=1:length(x)
                y = x(k);
                % Get precomputed values
                dxf = df(xfeats);
                
                %showNewton(x,x0,y,c,f,df,ddf);
                
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
        
        %         function showNewton(x, x0, ybase, c, f, df, ddf)
        %             % Fix an identical y for all iterations
        %             y = repmat(ybase,1,size(x,2));
        %             gs = optFun(x,y,x0,c,f,df,ddf);
        %             xs = x;
        %             for i=1:2
        %                 [g,dg] = optFun(x,y,x0,c,f,df,ddf);
        %                 xtmp = x;
        %                 x = x - g./dg;
        %                 plot(xs,gs,'r','LineWidth',2);
        %                 hold on;
        %                 plot([xtmp; x],[g; zeros(size(g))],x,0,'rx');
        %                 hold off;
        %                 title(sprintf('ybase=%f',ybase));
        %                 %axis([min(xs) max(xs) min(gs)*1.5 max(gs)*1.5]);
        %                 pause;
        %             end
        %         end
    end
end

