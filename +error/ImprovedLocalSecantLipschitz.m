classdef ImprovedLocalSecantLipschitz < error.BaseLocalLipschitzFunction
% ImprovedLocalSecantLipschitz: 
%
%
%
% @author Daniel Wirtz @date 2011-05-20
%
% @new{0,4,dw,2011-05-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
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
        
        % @propclass{optional}
        % @todo make docs, make dependent and re-generate if set to true when not used..
        UseCachedSecants = false;
    end

    properties(Access=private)
        oldxfeat = [];
        precomp = [];
    end
    
    methods
        
        function this = ImprovedLocalSecantLipschitz(bellfcn)
            this = this@error.BaseLocalLipschitzFunction(bellfcn);
            this.registerProps('NewtonTolerance','MaxNewtonIterations','UseCachedSecants');
        end
        
        function ci = evaluate(this, di, C, t, mu)%#ok
            b = this.bellfcn;
            x0 = b.x0;
            if isempty(this.oldxfeat) || any(isnan(this.oldxfeat)) || size(this.oldxfeat,2) ~= size(di,2)
                this.oldxfeat = x0+.2*sign(x0-di);
            end
            
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                if this.UseCachedSecants
                    xfeat = this.CachedModifiedNewton(di);
                else
                    xfeat = this.ModifiedNewton(this.oldxfeat,di);
                end
                ci = abs(b.evaluateD1(xfeat));
                this.oldxfeat = xfeat;
            else
                xfeat = ones(size(di))*x0;
                update = abs(di-x0) < C;
                % Choose suitable starting conditions if no old xfeat vectors
                % are available
                if any(update)
                    if this.UseCachedSecants
                        xfeat(update) = this.CachedModifiedNewton(di(update));
                    else
                        xfeat(update) = this.ModifiedNewton(this.oldxfeat(update),di(update));
                    end
                    this.oldxfeat(update) = xfeat(update);
                end
                left = di + C - xfeat < 0;
                right = di - C - xfeat > 0;
                center = ~left & ~right;
                % If C is too small we just take the derivative at this
                % spot
                if C < sqrt(eps)
                    ci(left | right) = abs(b.evaluateD1(di(left | right)));
                else
                    ci(left) = (b.evaluateScalar((di(left))) - b.evaluateScalar((di(left)+C))) / C;
                    ci(right) = (b.evaluateScalar((di(right)-C)) - b.evaluateScalar((di(right)))) / C;
                end
                ci(center) = abs(b.evaluateD1(xfeat(center)));
            end
        end
        
        function precompMaxSecants(this, maxDistance, num)
            if this.UseCachedSecants
                di = linspace(0,maxDistance,num);
                start = this.bellfcn.x0+.2*sign(this.bellfcn.x0-di);
                if KerMor.App.Verbose > 1
                    fprintf('Precomputing %d max local secants for BellFunction(Gamma=%f,x0=%f) in range [0, %f] ..\n',num,this.Gamma,this.x0,maxDistance);
                end
                xfeat = this.ModifiedNewton(start,di);
                t = general.collections.BinTree;
                if KerMor.App.Verbose > 1
                    fprintf('Building tree data structure..\n');
                end
                t.Insert(di,xfeat);
                this.precomp = t;
            end
        end
    end
    
    methods(Access=private)
        
        function xtmp = CachedModifiedNewton(this, y)
            if ~isempty(this.precomp)
                t = this.precomp;
                [l,u] = t.FindClosest(y);
                low = y < this.bellfcn.x0;
                xtmp = zeros(size(y));
                xtmp(low) = u(low);
                xtmp(~low) = l(~low);
            else
                error('When using cached secants the method precompMaxSecants has to be run first.');
            end
        end
        
        function xtmp = ModifiedNewton(this, xstart, y)
            xtmp = xstart+2*this.NewtonTolerance;
            x = xstart;
            
            b = this.bellfcn; x0 = b.x0; xR = b.xR;
            
            f = @(x)b.evaluateScalar(x);
            df = @(x)b.evaluateD1(x);
            ddf = @(x)b.evaluateD2(x);
            
            % Add eps at division as for zero nominator the expression is zero. (no error made)
            n0 = df(0) - (f(0)-f(y))./(eps-y); 
            nx0 = df(x0) - (f(x0)-f(y))./(x0-y+eps);
            nxr = df(xR) - (f(xR)-f(y))./(xR-y+eps);
            
            % Lim x->y n'(x) = \phi''(x) / 2 (De L'hospital)
            dn0 = ddf(0) - n0./-y;
            dn0(isnan(dn0)) = ddf(0)/2;
            dnx0 = ddf(x0) - nx0./(x0-y);
            dnx0(isinf(dnx0)) = 0; % == ddf(x0)/2;
            dnxr = ddf(xR) - nxr./(xR-y);
            dnxr(isinf(dnxr)) = ddf(xR)/2;
            
            cnt = 0;
            
            %conv = 1:length(x);
            %finished = zeros(size(x));
            %while ~isempty(conv) && cnt < b.MaxNewtonIterations
            while any(abs(xtmp-x) > this.NewtonTolerance) && cnt < this.MaxNewtonIterations
                
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
                
                g = zeros(size(x));
                dg = g;
                
                a = y < x0;
                pi = a & x < x0 | ~a & x > x0;
                pr = a & x > xR;
                pl = ~a & x < 0;
                std = ~(pi | pl | pr);

                % Standard case
                xs = x(std); ys = y(std); 
                g(std) = df(xs) - (f(xs)-f(ys))./(xs-ys);
                g(isnan(g)) = 0;
                dg(std) = ddf(xs) - g(std)./(xs-ys);
                ninf = isinf(dg);
                dg(ninf) = ddf(xs(ninf))/2;
                
                % Penalty case
                p = ~std;
                d = zeros(size(x)); c = d;
                c(pi) = -dnx0(pi).^2 ./ (4*nx0(pi));
                d(pi) = 2*nx0(pi)./dnx0(pi) - x0;
                c(pl) = -dn0(pl).^2 ./ (4*n0(pl));
                d(pl) = 2*n0(pl)./dn0(pl);
                c(pr) = dnxr(pr).^2 ./ (4*nxr(pr));
                d(pr) = 2*nxr(pr)./dnxr(pr) - xR;
                
                g(p) = c(p) .* (x(p) + d(p)).^2;
                dg(p) = 2 * c(p) .* x(p);
                
                xtmp = x;
                x = x - g./dg;
                cnt = cnt + 1;
                
                %if cnt > this.MaxNewtonIterations*.5
                    %showNewton;
                %end
            end
            if cnt == this.MaxNewtonIterations
                error('Bellfunction->ModifiedNewton: Max iterations of %d reached',this.MaxNewtonIterations);
            end
            
            function showNewton %#ok<*DEFNU>
                % Fix an identical y for all iterations
                si = 200;
                for i=1:5
                    
                    X = linspace(.5*min(min([x; xtmp],[],1)),2*max(max([x; xtmp],[],1)),si);
                    gs = zeros(length(y),si);
                    for idx = 1:length(y)
                        gs(idx,:) = error.ImprovedLocalSecantLipschitz.optFun(X,repmat(y(idx),size(X)),f,df,ddf,this.x0);
                    end
                    
                    [g,dg] = error.ImprovedLocalSecantLipschitz.optFun(x,y,f,df,ddf,this.x0);
                    xtmp = x;
                    x = x - g./dg;
                    %plot(xs,gs,'r');%,'LineWidth',2
                    figure(1);
                    subplot(1,2,1);
                    plot(X,gs,'b',X,0,'black-');
                    hold on;
                    plot([xtmp; x],[g; zeros(size(g))],x,0,'rx',xtmp,g,'rx');
                    hold off;
                    axis tight;
                    
                    subplot(1,2,2);
                    X = linspace(0-this.xR,2*this.xR,si);
                    gs = zeros(length(y),si);
                    for idx = 1:length(y)
                        gs(idx,:) = kernels.BellFunction.optFun(X,repmat(y(idx),size(X)),f,df,ddf,this.x0,this.PenaltyFactor);
                    end
                    plot(X,gs);
                    hold on;
                    plot(X,0,'black-');
                    hold off;
                    pause;
                end
            end
        end
    end
    
    methods(Static)
        function xFeatDemo(Gamma)
            if nargin == 0
                Gamma = 2;
            end
            
            bell = kernels.GaussKernel(Gamma);
            lfun = error.ImprovedLocalSecantLipschitz(bell);
            lfun.NewtonTolerance = 1e-3;
            h = figure(1);
            
            maxX = Gamma*3;
            x = linspace(0,maxX,7000);
            x = union(x,[bell.x0,bell.xR]);
            
            %xstart = x0*(1+.5*sign(x0-x));
            xstart = repmat(bell.x0/2,1,size(x,2));
            
            xfeats = lfun.ModifiedNewton(xstart,x);
            
            f = @bell.evaluateScalar;
            df = @bell.evaluateD1;
            ddf = @bell.evaluateD2;
            
            %kernels.BellFunction.showNewton(x,bell.x0,y,f,df,ddf);

            for k=1:length(x)
                y = x(k);

                % Get precomputed values
                dxf = df(xfeats);
                
                %kernels.BellFunction.showNewton(x,bell.x0,y,f,df,ddf);
                
                % Get optimization function
                [opt,dummy,inner,l0,gxr] = error.ImprovedLocalSecantLipschitz.optFun(x,repmat(y,1,size(x,2)),f,df,ddf,bell.x0);
                std = ~inner & ~l0 & ~gxr;
                
                % Bound plot
                hlp = max(max(abs(df(x))),max(abs(f(x))));
                %hlp = (maxX-y)*dxf(k); % also used in plot of tangent!
                opt(opt < -.2*hlp) = -.2*hlp;
                %hlp2 = max(1.5,max(df(x)));
                opt(opt > hlp) = hlp;
                
                % Plot!
                plot(x,f(x),'r', x,abs(df(x)),'m', x,abs(dxf),'g',[0 y maxX],f(y) + [-y*dxf(k) 0 hlp],...
                    'g--',bell.x0,f(bell.x0),'black.',y,f(y),'r.',bell.xR,f(bell.xR),'b.');
                hold on;
                st = {};
                if any(std)
                    plot(x(std),opt(std),'b');
                    st{end+1} = 'Opt. target'; %#ok<*AGROW>
                end
                if any(inner)
                    plot(x(inner),opt(inner),'b--');
                    st{end+1} = 'penalty inner';
                end
                if any(l0)
                    plot(x(l0),opt(l0),'b-.');
                    st{end+1} = 'penalty leq 0';
                end
                if (any(gxr))
                    plot(x(gxr),opt(gxr),'b-.');
                    st{end+1} = 'penalty geq xr';
                end
                legend('Gaussian','Absolute Gaussian derivative','Max local derivatives',...
                    'Max local gradient around x0','x_0','y','x_R',st{:});
                plot([xfeats(k) xfeats(k) y],[f(xfeats(k)) abs(dxf(k)) abs(dxf(k))],'blackx',x,0,'black');
                hold off;
                pause;
            end
            
            close(h);
            
        end
        
        function [g,dg,p1,p2,p4] = optFun(x,y,f,df,ddf,x0)
            g = zeros(size(x));
            dg = g;
            
            xr = f(0)*x0 / (f(0)-f(x0));
            n0 = df(0) - (f(0)-f(y))./(-y+eps);
            nx0 = df(x0) - (f(x0)-f(y))./(x0-y+eps);
            nxr = df(xr) - (f(xr)-f(y))./(xr-y+eps);
            
            dn0 = ddf(0) - n0./-y;
            dn0(isnan(dn0)) = ddf(0)/2;
            dnx0 = ddf(x0) - nx0./(x0-y);
            dnx0(isinf(dnx0)) = 0; % == ddf(this.x0)/2;
            dnxr = ddf(xr) - nxr./(xr-y);
            dnxr(isinf(dnxr)) = ddf(xr)/2;
            
            a = y <= x0;
%             std = true(1,length(x));
%             p1 = ~std; p2 = p1; p4 = p1;
            std = (a & x < xr & x > x0) | (~a & 0 < x & x < x0);
            p1 = a & x <= x0 | ~a & x > x0;
            p4 = a & x >= xr;
            p2 = ~a & x <= 0;
            
            % Standard case
            xs = x(std); ys = y(std);
            g(std) = df(xs) - (f(xs)-f(ys))./(xs-ys);
            if any(isnan(g))
                keyboard;
            end
            g(isnan(g)) = 0;
            dg(std) = ddf(xs) - g(std)./(xs-ys);
            iinf = isinf(dg);
            dg(iinf) = ddf(xs(iinf))/2;
            
            if any(p1)
%                 g(p1) = nx0(p1) - p*(x(p1)-x0).^2;
%                 dg(p1) = -2*p*(x(p1)-x0);
                c = dnx0(p1).^2 ./ (4*nx0(p1)); %sign(dnx0(p1)).*
                g(p1) = c.*(x(p1) + 2*nx0(p1)./dnx0(p1) - x0).^2;
                dg(p1) = 2*c.*x(p1);
%                 xb = (xr+x0)/2;
%                 g(p1) = nx0(p1) - dnx0(p1).*(xb/2 + (x(p1)+xb).^2 / (2*xb));
%                 dg(p1) = dnx0(p1).*x(p1) + 1;
            end
            if any(p2) % FIX!
%                 g(p2) = n0(p2) + p*x(p2).^2;
%                 dg(p2) = 2*p*x(p2);
%                 xb = -xr;
%                 g(p2) = n0(p2) - dn0(p2).*(xb/2 + (x(p2)+xb).^2 / (2*xb));
%                 dg(p2) = dn0(p2).*x(p2) + 1;
                c = dn0(p2).^2 ./ (4*n0(p2)); % sign(dn0(p2)).*
                g(p2) = c.*(x(p2) + 2*n0(p2)./dn0(p2) - 0).^2;
                dg(p2) = 2*c.*x(p2);
            end
            if any(p4)
%                 g(p4) = nxr(p4) + p*(x(p4)-xr).^2;
%                 dg(p4) = 2*p*(x(p4)-xr);
%                 xb = 1.5*xr;
%                 g(p4) = nxr(p4) - dnxr(p4).*(xb/2 + (x(p4)+xb).^2 / (2*xb));
%                 dg(p4) = dnxr(p4).*x(p4) + 1;
                c = sign(dnxr(p4)).*dnxr(p4).^2 ./ (4*nxr(p4));
                g(p4) = c.*(x(p4) + 2*nxr(p4)./dnxr(p4) - xr).^2;
                dg(p4) = 2*c.*x(p4);
            end
        end
        
        function showNewton(x, x0, ybase, f, df, ddf)
            % Fix an identical y for all iterations
            y = repmat(ybase,1,size(x,2));
            gs = kernels.BellFunction.optFun(x,y,f,df,ddf,x0);
            xs = x;
            for i=1:5
                [g,dg] = error.ImprovedLocalSecantLipschitz.optFun(x,y,f,df,ddf,x0);
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