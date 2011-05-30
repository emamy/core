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
        oldrs = [];
        precomp = [];
    end
    
    methods
        
        function this = ImprovedLocalSecantLipschitz(bellfcn)
            this = this@error.BaseLocalLipschitzFunction(bellfcn);
            this.registerProps('NewtonTolerance','MaxNewtonIterations','UseCachedSecants');
        end
        
        function ci = evaluate(this, di, C, t, mu)%#ok
            b = this.bellfcn;
            r0 = b.r0;
            if isempty(this.oldrs) || any(isnan(this.oldrs)) || size(this.oldrs,2) ~= size(di,2)
                this.oldrs = r0+.2*sign(r0-di);
            end
            
            % Consider C=Inf case separately for speed reasons
            if isinf(C)
                if this.UseCachedSecants
                    rs = this.CachedModifiedNewton(di);
                else
                    rs = this.ModifiedNewton(this.oldrs,di);
                end
                ci = abs(b.evaluateD1(rs));
                this.oldrs = rs;
            else
                rs = ones(size(di))*r0;
                update = abs(di-r0) < C;
                % Choose suitable starting conditions if no old rs vectors
                % are available
                if any(update)
                    if this.UseCachedSecants
                        rs(update) = this.CachedModifiedNewton(di(update));
                    else
                        rs(update) = this.ModifiedNewton(this.oldrs(update),di(update));
                    end
                    this.oldrs(update) = rs(update);
                end
                left = di + C - rs < 0;
                right = di - C - rs > 0;
                center = ~left & ~right;
                % If C is too small we just take the derivative at this
                % spot
                if C < sqrt(eps)
                    ci(left | right) = abs(b.evaluateD1(di(left | right)));
                else
                    ci(left) = (b.evaluateScalar((di(left))) - b.evaluateScalar((di(left)+C))) / C;
                    ci(right) = (b.evaluateScalar((di(right)-C)) - b.evaluateScalar((di(right)))) / C;
                end
                ci(center) = abs(b.evaluateD1(rs(center)));
            end
        end
        
        function precompMaxSecants(this, maxDistance, num)
            if this.UseCachedSecants
                di = linspace(0,maxDistance,num);
                start = this.bellfcn.r0+.2*sign(this.bellfcn.r0-di);
                if KerMor.App.Verbose > 1
                    fprintf('Precomputing %d max local secants for BellFunction(Gamma=%f,r0=%f) in range [0, %f] ..\n',num,this.Gamma,this.x0,maxDistance);
                end
                rs = this.ModifiedNewton(start,di);
                t = general.collections.BinTree;
                if KerMor.App.Verbose > 1
                    fprintf('Building tree data structure..\n');
                end
                t.Insert(di,rs);
                this.precomp = t;
            end
        end
    end
    
    methods(Access=private)
        
        function rs = CachedModifiedNewton(this, s)
            if ~isempty(this.precomp)
                t = this.precomp;
                [l,u] = t.FindClosest(s);
                low = s < this.bellfcn.r0;
                rs = zeros(size(y));
                rs(low) = u(low);
                rs(~low) = l(~low);
            else
                error('When using cached secants the method precompMaxSecants has to be run first.');
            end
        end
        
        function rtmp = ModifiedNewton(this, rstart, s)
            rtmp = rstart+2*this.NewtonTolerance;
            r = rstart;
            
            b = this.bellfcn; r0 = b.r0; rR = b.rR;
            
            % Add eps at division as for zero nominator the expression is zero. (no error made)
            n(1,:) = b.evaluateD1(0) - (b.evaluateScalar(0)-b.evaluateScalar(s))./(eps-s); 
            n(2,:) = b.evaluateD1(r0) - (b.evaluateScalar(r0)-b.evaluateScalar(s))./(r0-s+eps);
            n(3,:) = b.evaluateD1(rR) - (b.evaluateScalar(rR)-b.evaluateScalar(s))./(rR-s+eps);
            
            % Lim x->y n'(x) = \phi''(x) / 2 (De L'hospital)
            dn(1,:) = b.evaluateD2(0) - n(1,:)./-s;
            dn(1,isnan(dn(1,:))) = b.evaluateD2(0)/2;
            dn(2,:) = b.evaluateD2(r0) - n(2,:)./(r0-s);
            dn(2,isinf(dn(2,:))) = 0; % == ddf(x0)/2;
            dn(3,:) = b.evaluateD2(rR) - n(3,:)./(rR-s);
            dn(3,isinf(dn(3,:))) = b.evaluateD2(rR)/2;
            
            cnt = 0;
            
            %conv = 1:length(x);
            %finished = zeros(size(x));
            %while ~isempty(conv) && cnt < b.MaxNewtonIterations
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
                
                [g,dg] = this.optFun(r,s,n,dn,b);
                
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
%                         gs(idx,:) = error.ImprovedLocalSecantLipschitz.optFun(X,repmat(y(idx),size(X)),f,df,ddf,this.x0);
%                     end
%                     
%                     [g,dg] = error.ImprovedLocalSecantLipschitz.optFun(x,y,f,df,ddf,this.x0);
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
        
        function [g,dg,pi,pl,pr] = optFun(this,r,s,n,dn,b)%#ok
            g = zeros(size(r));
            d = g; c = d; dg = g;
            
            a = s < b.r0;
            pi = a & r < b.r0 | ~a & r > b.r0;
            pr = a & r > b.rR;
            pl = ~a & r < 0;
            std = ~(pi | pl | pr);
            
            % Standard case
            rs = r(std); ss = s(std);
            g(std) = b.evaluateD1(rs) - (b.evaluateScalar(rs)-b.evaluateScalar(ss))./(rs-ss);
            g(isnan(g)) = 0;
            dg(std) = b.evaluateD2(rs) - g(std)./(rs-ss);
            ninf = isinf(dg);
            dg(ninf) = b.evaluateD2(rs(ninf))/2;
            
            % Penalty case
            p = ~std;
            c(pi) = dn(2,pi).^2 ./ (4*n(2,pi));
            d(pi) = 2*n(2,pi)./dn(2,pi) - b.r0;
            c(pl) = dn(1,pl).^2 ./ (4*n(1,pl));
            d(pl) = 2*n(1,pl)./dn(1,pl);
            c(pr) = dn(3,pr).^2 ./ (4*n(3,pr));
            d(pr) = 2*n(3,pr)./dn(3,pr) - b.rR;
            
            g(p) = c(p) .* (r(p) + d(p)).^2;
            dg(p) = 2 * c(p) .* (r(p) + d(p));
            % Avoid exact matches!
            g(dg == 0) = 0;
            dg(dg == 0) = 1;
        end
    end
    
    methods(Static)
        function rsDemo(Gamma)
            if nargin == 0
                Gamma = 2;
            end
            
            b = kernels.GaussKernel(Gamma);
            r0 = b.r0; rR = b.rR;
            lfun = error.ImprovedLocalSecantLipschitz(b);
            lfun.NewtonTolerance = 1e-3;
            
            h = figure(1);
            plotrs = false;
            lw = 2;
            
            maxR = Gamma*3;
            r = linspace(0,maxR,70);
            r = union(r,[r0,rR]);
            
            rstart = repmat(r0/2,1,size(r,2));
            rss = lfun.ModifiedNewton(rstart,r);
            
            f = @b.evaluateScalar;
            df = @b.evaluateD1;
            
            %kernels.BellFunction.showNewton(x,bell.x0,y,f,df,ddf);
            m = length(r);
            for k=1:m
                s = r(k);

                % Get precomputed values
                drf = abs(df(rss));
                
                %kernels.BellFunction.showNewton(x,bell.x0,y,f,df,ddf);
                
                % Add eps at division as for zero nominator the expression is zero. (no error made)
                n(1,1:m) = b.evaluateD1(0) - (b.evaluateScalar(0)-b.evaluateScalar(s))./(eps-s); 
                n(2,1:m) = b.evaluateD1(r0) - (b.evaluateScalar(r0)-b.evaluateScalar(s))./(b.r0-s+eps);
                n(3,1:m) = b.evaluateD1(rR) - (b.evaluateScalar(rR)-b.evaluateScalar(s))./(b.rR-s+eps);

                % Lim x->y n'(x) = \phi''(x) / 2 (De L'hospital)
                dn(1,1:m) = b.evaluateD2(0) - n(1,:)./-s;
                dn(1,isnan(dn(1,:))) = b.evaluateD2(0)/2;
                dn(2,1:m) = b.evaluateD2(r0) - n(2,:)./(r0-s);
                dn(2,isinf(dn(2,:))) = 0; % == ddf(x0)/2;
                dn(3,1:m) = b.evaluateD2(rR) - n(3,:)./(rR-s);
                dn(3,isinf(dn(3,:))) = b.evaluateD2(rR)/2;
                
                % Get optimization function
                [opt,dummy,inner,l0,gxr] = lfun.optFun(r,repmat(s,1,size(r,2)),n,dn,b);
                std = ~inner & ~l0 & ~gxr;
                                
                % f and df
                plot(r,f(r),'r',r,abs(df(r)),'m','LineWidth',lw);
                st={};
                st(1:2) = {'\phi', '|\phi''|'}; %#ok<*AGROW>
                hold on;
                
                % Plot abs(df(rs)) for all s
                if plotrs
                    plot(r,abs(drf),'g','LineWidth',lw);
                    st{end+1} = '|\phi''(r_s)| = |S_s(r_s))| \forall s';
                end
                % Plot f tangent for max secant gradient
                plot([0 maxR],f(s) + drf(k)*[s (s-maxR)],'g--','LineWidth',lw);
                st{end+1} = '\phi(r_s) + \phi''(r_s)(r_s-s)';
                % Plot 
                plot(r0,f(r0),'blacks','MarkerSize',6,'MarkerFaceColor','black');
                plot(rR,f(rR),'bs','MarkerSize',6,'MarkerFaceColor','blue');
                plot(s,f(s),'ro','MarkerSize',6,'MarkerFaceColor','red');
                plot(rss(k),f(rss(k)),'mo','MarkerSize',6,'MarkerFaceColor','magenta');
                st(end+1:end+4) = {'\phi(r_0)','\phi(r_R)','\phi(s)','\phi(r_s)'};
                if any(std)
                    plot(r(std),opt(std),'b','LineWidth',lw);
                    st{end+1} = 'n(r)'; 
                end
                if any(inner)
                    plot(r(inner),opt(inner),'b--','LineWidth',lw);
                    st{end+1} = 'c(r_0)(r+d(r_0))^2';
                end
                if any(l0)
                    plot(r(l0),opt(l0),'b-.','LineWidth',lw);
                    st{end+1} = 'c(0)(r+d(0))^2';
                end
                if (any(gxr))
                    plot(r(gxr),opt(gxr),'b-.','LineWidth',lw);
                    st{end+1} = 'c(r_R)(r+d(r_R))^2';
                end
                plot([rss(k) rss(k)],[f(rss(k)) abs(drf(k))],'--black','LineWidth',lw);
                if plotrs
                    plot([rss(k) s],[abs(drf(k)) abs(drf(k))],'--black','LineWidth',lw);
                end
                st{end+1} = '\phi(r_s) <-> d\phi(rs) = |S_s(r_s)|';
                % Plot line through zero
                plot(r,0,'black');
                axis([0 maxR -.5 f(r0)+abs(df(r0))*r0]);
                hold off;
                
                legend(st{:});
                pause;
            end
            
            close(h);
            
        end
        
%         function showNewton(x, x0, ybase, f, df, ddf)
%             % Fix an identical y for all iterations
%             y = repmat(ybase,1,size(x,2));
%             gs = kernels.BellFunction.optFun(x,y,f,df,ddf,x0);
%             xs = x;
%             for i=1:5
%                 [g,dg] = error.ImprovedLocalSecantLipschitz.optFun(x,y,f,df,ddf,x0);
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
    
end