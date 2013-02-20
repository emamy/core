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
            p = ~std;
            c(pi) = dn(2,pi).^2 ./ (4*n(2,pi));
            d(pi) = 2*n(2,pi)./dn(2,pi) - this.r0;
            c(pl) = dn(1,pl).^2 ./ (4*n(1,pl));
            d(pl) = 2*n(1,pl)./dn(1,pl);
            c(pr) = dn(3,pr).^2 ./ (4*n(3,pr));
            d(pr) = 2*n(3,pr)./dn(3,pr) - this.rm;
            
            g(p) = c(p) .* (r(p) + d(p)).^2;
            dg(p) = 2 * c(p) .* (r(p) + d(p));
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
        function KMdemo_rsDemoImages(b)
            % Produces the demo images for rs positions
            if nargin == 0
                b = kernels.GaussKernel(2);
            end
            
            r0 = b.r0; rm = b.rm; %#ok<*PROP>
            
            f = @b.evaluateScalar;
            df = @b.evaluateD1;
            
            
            lw = 1;
            fs = 16;
            
            maxR = rm*2;
            r = linspace(0,maxR,70);
            
            %alls = [r0 rm+r0]/2;
            alls = [.2*r0 1.5*rm];
            rs = b.ModifiedNewton([r0 r0],alls);
            drs = abs(df(rs));
            
            for k=1:2
                h = figure(1);
                s = alls(k);
                % f and df
                plot(r,f(r),'r','LineWidth',lw); %,r,abs(df(r)),'m'
                hold on;
                
                % Plot f tangent for max secant gradient
                plot([0 maxR],f(s) + drs(k)*[s (s-maxR)],'b--','LineWidth',lw);
                % Plot & text
                xoff = maxR*.001;
                yoff = .05;
                plot([r0 r0],[0 f(r0)],'blacks','MarkerSize',6,'MarkerFaceColor','black');
                text(r0,-yoff,'  r_0','FontSize',fs);
                text(r0+xoff,f(r0)+yoff,'  \phi(r_0)','FontSize',fs);
                
                plot([rm rm],[0 f(rm)],'bs','MarkerSize',6,'MarkerFaceColor','blue');
                text(rm,-yoff,'  r_m','FontSize',fs);
                text(rm+xoff,f(rm)+yoff,'  \phi(r_m)','FontSize',fs);
                
                plot([s s],[0 f(s)],'ro','MarkerSize',6,'MarkerFaceColor','red');
                text(s,-yoff,'  s','FontSize',fs);
                text(s+xoff,f(s)+yoff,'  \phi(s)','FontSize',fs);
                
                plot([rs(k) rs(k)],[0 f(rs(k))],'mo','MarkerSize',6,'MarkerFaceColor','magenta');
                text(rs(k),-yoff,'  r_s','FontSize',fs);
                text(rs(k)+xoff,f(rs(k))+yoff,'  \phi(r_s)','FontSize',fs);
                
                plot(r,0,'black');
               
                axis([0 maxR -.5 f(r0)+abs(df(r0))*r0]);
                hold off;
                
                lh = legend('\phi','\phi(r_s) + \phi''(r_s)(x-r_s)');%, '|\phi''|'
                set(lh,'FontSize',fs);
                title(sprintf('Bell function demo, s=%f, rs=%f',s,rs(k)),'FontSize',fs);
                pause;
                Utils.saveFigure(h,sprintf('rsDemo%d',k),'eps');
                %saveas(h,sprintf('rsDemo%ddirect.eps',k),'eps2c');
                close(h);
            end
            
        end
        
        function penaltyDemo(Gamma)
            if nargin == 0
                Gamma = 2;
            end
            
            b = kernels.GaussKernel(Gamma);
            r0 = b.r0; rm = b.rm;
            lfun = error.lipfun.ImprovedLocalSecantLipschitz(b);
            b.NewtonTolerance = 1e-3;
            lfun.prepareConstants;
            
            h = figure(1);
            plotrs = false;
            lw = 1;
            
            maxR = Gamma*3;
            r = linspace(0,maxR,70);
            r = union(r,[r0,rm]);
            
            rstart = repmat(r0/2,1,size(r,2));
            rss = b.ModifiedNewton(rstart,r);
            
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
                n(3,1:m) = b.evaluateD1(rm) - (b.evaluateScalar(rm)-b.evaluateScalar(s))./(b.rm-s+eps);

                % Lim x->y n'(x) = \phi''(x) / 2 (De L'hospital)
                dn(1,1:m) = b.evaluateD2(0) - n(1,:)./-s;
                dn(1,isnan(dn(1,:))) = b.evaluateD2(0)/2;
                dn(2,1:m) = b.evaluateD2(r0) - n(2,:)./(r0-s);
                dn(2,isinf(dn(2,:))) = 0; % == ddf(x0)/2;
                dn(3,1:m) = b.evaluateD2(rm) - n(3,:)./(rm-s);
                dn(3,isinf(dn(3,:))) = b.evaluateD2(rm)/2;
                
                % Get optimization function
                [opt,~,inner,l0,gxr] = b.optFun(r,repmat(s,1,size(r,2)),n,dn);
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
                plot([0 maxR],f(s) + drf(k)*[s (s-maxR)],'k--','LineWidth',lw);
                st{end+1} = '\phi(r_s) + \phi''(r_s)(x-r_s)';
                % Plot 
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
                    st{end+1} = 'c(r_m)(r+d(r_m))^2';
                end
                plot([rss(k) rss(k)],[f(rss(k)) abs(drf(k))],'--black','LineWidth',lw);
                if plotrs
                    plot([rss(k) s],[abs(drf(k)) abs(drf(k))],'--black','LineWidth',lw);
                end
                st{end+1} = '\phi(r_s) <-> d\phi(rs) = |S_s(r_s)|';
                
                fs = 16;
                off = .04; off2 = .07;
                plot([r0 r0],[0 f(r0)],'blacks','MarkerSize',6,'MarkerFaceColor','black');
                text(r0+off,f(r0)+off,'\phi(r_0)','FontSize',fs);
                text(r0,-off2,'r_0','FontSize',fs);
                plot([rm rm],[0 f(rm)],'bs','MarkerSize',6,'MarkerFaceColor','blue');
                text(rm+off,f(rm)+off,'\phi(r_m)','FontSize',fs);
                text(rm,-off2,'r_m','FontSize',fs);
                plot([s s],[0 f(s)],'ro','MarkerSize',6,'MarkerFaceColor','red');
                text(s+off,f(s)+off,'\phi(s)','FontSize',fs);
                text(s,-off2,'s','FontSize',fs);
                plot([rss(k) rss(k)],[0 f(rss(k))],'mo','MarkerSize',6,'MarkerFaceColor','magenta');
                text(rss(k)+off,f(rss(k))+off,'\phi(r_s)','FontSize',fs);
                text(rss(k),-off2,'r_s','FontSize',fs);
                %st(end+1:end+4) = {'\phi(r_0)','\phi(r_m)','\phi(s)','\phi(r_s)'};
                
                
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

