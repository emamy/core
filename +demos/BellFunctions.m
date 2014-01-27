classdef BellFunctions
    % BellFunctions: Demos regarding the Bell function local Lipschitz
    % estimations
    %
    % We kindly refer to the article \cite WH12 for introduction of bell
    % functions and their theory.
    %
    % @author Daniel Wirtz @date 2013-08-23
    %
    % @new{0,7,dw,2013-08-23} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    methods(Static)
        function pm = SecantGradientPlots(bfun)
            % Produces the demo images for maximum secant gradient
            % positions.
            %
            % For details on the implemented methodology please see
            % @cite WH12
            %
            % Returns a PlotManager handle if the plots are to be saved.
            %
            % Parameters:
            % bfun: A bell function @type kernels.BellFunction @default
            % kernels.GaussKernel(2)
            %
            % Return values:
            % pm: A PlotManager instance @type PlotManager
            if nargin == 0
                bfun = kernels.GaussKernel(2);
            end
            
            r0 = bfun.r0; rm = bfun.rm; %#ok<*PROP>
            
            f = @bfun.evaluateScalar;
            df = @bfun.evaluateD1;
            
            lw = 1;
            fs = 16;
            
            maxR = rm*2;
            r = linspace(0,maxR,70);
            
            %alls = [r0 rm+r0]/2;
            alls = [.2*r0 1.5*rm];
            rs = bfun.ModifiedNewton([r0 r0],alls);
            drs = abs(df(rs));
            
            pm = PlotManager;
            pm.LeaveOpen = true;
            pm.SaveFormats = {'eps'};
            
            for k=1:2
                s = alls(k);
                h = pm.nextPlot(sprintf('secant_demo_%d',k),...
                    sprintf('Maximum secant gradient selection for bell functions, s=%f, rs=%f',s,rs(k)));

                % f and df
                plot(h,r,f(r),'r','LineWidth',lw); %,r,abs(df(r)),'m'
                hold(h, 'on');
                
                % Plot f tangent for max secant gradient
                plot(h,[0 maxR],f(s) + drs(k)*[s (s-maxR)],'b--','LineWidth',lw);
                % Plot & text
                xoff = maxR*.001;
                yoff = .05;
                plot(h,[r0 r0],[0 f(r0)],'blacks','MarkerSize',6,'MarkerFaceColor','black');
                text(r0,-yoff,'  r_0','FontSize',fs,'Parent',h);
                text(r0+xoff,f(r0)+yoff,'  \phi(r_0)','FontSize',fs,'Parent',h);
                
                plot(h,[rm rm],[0 f(rm)],'bs','MarkerSize',6,'MarkerFaceColor','blue');
                text(rm,-yoff,'  r_m','FontSize',fs,'Parent',h);
                text(rm+xoff,f(rm)+yoff,'  \phi(r_m)','FontSize',fs,'Parent',h);
                
                plot(h,[s s],[0 f(s)],'ro','MarkerSize',6,'MarkerFaceColor','red');
                text(s,-yoff,'  s','FontSize',fs,'Parent',h);
                text(s+xoff,f(s)+yoff,'  \phi(s)','FontSize',fs,'Parent',h);
                
                plot(h,[rs(k) rs(k)],[0 f(rs(k))],'mo','MarkerSize',6,'MarkerFaceColor','magenta');
                text(rs(k),-yoff,'  r_s','FontSize',fs,'Parent',h);
                text(rs(k)+xoff,f(rs(k))+yoff,'  \phi(r_s)','FontSize',fs,'Parent',h);
                
                plot(h,r,0,'black');
                
                axis(h,[0 maxR -.5 f(r0)+abs(df(r0))*r0]);
                
                lh = legend(h,'\phi','\phi(r_s) + \phi''(r_s)(x-r_s)');%, '|\phi''|'
                set(lh,'FontSize',fs);
            end
        end
        
        function NewtonPenalty(Gamma)
            % Demonstrates the error estimator penalized newton function
            % and maximum secant gradients
            %
            % Parameters:
            % Gamma: The `\gamma` value to use for the Gaussian bell
            % function. @type double @default 2
            if nargin < 1
                Gamma = 2;
            end
            
            b = kernels.GaussKernel(Gamma);
            r0 = b.r0; rm = b.rm;
            lfun = error.lipfun.ImprovedLocalSecantLipschitz(b);
            b.NewtonTolerance = 1e-3;
            lfun.prepareConstants;
            
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
            fh = figure;
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
                
                title(sprintf('Maximum secant gradient r_s=%g for s=%g',rss(k),s));
                
                axis([0 maxR -.5 f(r0)+abs(df(r0))*r0]);
                hold off;
                
                legend(st{:});
                pause(1);
                % Stop if figure has been closed
                if ~ishandle(fh)
                    return;
                end
            end
        end
        
        function LocalLipschitzDemo( x0, C )
            %ERRORESTDEMO Demo for the monotone radial basis functions error estimator.
            
            h = figure(1);
            dt = 0.1;
            b = 1;
            
            if nargin < 2
                Cint = .1:.5:3;
                %Cint = 2;
                if nargin < 1
                    x0 = sqrt(b/2)+2;
                end
            else
                Cint = C;
            end
            Cint(end+1) = Inf;
            
            x = 0:dt:8.5;
            x = union(x,x0);
            
            bfun = kernels.GaussKernel(b);
            
            f = @(x)bfun.evaluateScalar(x);
            df = @(x)bfun.evaluateD1(x);
            ddf = @(x)bfun.evaluateD2(x);
            
            fx = f(x);
            maxR = bfun.r0;
            
            % precompute xfeats-vector
            %xfeats = newton(maxR+sign(maxR-x),x,f,df,ddf,...
            %    K.NewtonTolerance,K.r0,K.PenaltyFactor);
            
            for C = Cint
                maxder = ones(size(x))*abs(df(bfun.r0));
                LGL = maxder;
                LSE = maxder;
                minder = zeros(size(x));
                maxd = zeros(size(x));
                for i = 1:length(x)
                    d = abs(x(i));
                    
                    % Less efficient for non-secant case, but shows
                    % advantage of estimation on [maxR, xfeat] more nicely
                    if d - C - maxR > 0
                        LGL(i) = abs(df(d-C));
                        LSE(i) = (f(d-C) - f(d)) / C;
                    elseif d + C - maxR < 0
                        LGL(i) = abs(df(d+C));
                        LSE(i) = (f(d) - f(d+C)) / C;
                    else
                        xfeat = bfun.ModifiedNewton(bfun.r0/2,d);
                        if d + C - xfeat < 0
                            LSE(i) = (f(d) - f(d+C)) / C;
                        elseif d - C - xfeat > 0
                            LSE(i) = (f(d-C) - f(d)) / C;
                        else
                            LSE(i) = abs(df(xfeat));
                        end
                    end
                    
                    
                    % More efficient method (better GLE)
%                     xfeat = maxR;
%                     if abs(d-maxR) < C
%                         xfeat = bfun.ModifiedNewton(bfun.r0/2,d);
%                     end
%                     if d + C - xfeat < 0
%                         LGL(i) = abs(df(d+C));
%                         LSE(i) = (f(d) - f(d+C)) / C;
%                     elseif d - C - xfeat > 0
%                         LGL(i) = abs(df(d-C));
%                         LSE(i) = (f(d-C) - f(d)) / C;
%                     else
%                         LSE(i) = abs(df(xfeat));
%                     end
                    
                    
                    % Minimum possible derivative
                    di = abs(x - x(i));
                    idx = di <= C;
                    minder(i) = max(abs(fx(idx) - exp(-x(i).^2/b)) ./ di(idx));
                    
                    % maxD-Comp
                    d = (x(i)-x0);
                    maxd(i) = (fx(i)-exp(-x0.^2/b))/d;
                    
                end
                
                ph = plot(x,f(x),'r',x,abs(df(x)),'r--',x,LGL,'b',x,LSE,'g');
                set(ph,'LineWidth',2);
                hold on;
                title(sprintf('C=%f',C));
                
                skip = round(1/(3*dt));
                plot(x(1:skip:end),minder(1:skip:end),'m^','MarkerSize',6);
                legend('Gaussian','Abs. Gaussian derivative',...
                    'Local Gradient estimate (LGL)','Local secant estimate (LSL)',...
                    'Minimal Lipschitz constant');
                
                % C-range plot
                plot([x0-C x0-C+eps], [0 1], 'black--',[x0 x0+eps], [0 1], 'black',[x0+C x0+C+eps], [0 1], 'black--');
                                
                % maxR-Stelle (sqrt(b/2))
                plot(x,0,'black');
                plot(maxR,bfun.evaluateScalar(maxR),'r.','MarkerSize',15);
                
                % maxderivative at x0-Plot
%                 plot(x,maxd,'black');
%                 plot(x,max(maxd),'black');
                
                %     x0idx = find(x==x0,1);
                %     stx0 = -e3(x0idx);
                %     stxfeat = df(xfeats(x0idx));
                %     x0feat = xfeats(x0idx);
                % xfeats-gerade
                %     plot([0 x0 5],[-x0*stxfeat+f(x0) f(x0) (5-x0)*stxfeat+f(x0)],'g--',x0,abs(stxfeat),'green.');
                %     plot([0 x0 5],[-x0*stx0+f(x0) f(x0) (5-x0)*stx0+f(x0)],'magenta--',x0,abs(stxfeat),'magenta.');
                %     plot([x0feat x0feat+eps],[0 1],'y');
                %plot([x0 xfeats(x0idx)],[f(x0) sign(maxR-x0)*(x0-x0feat)*stx0+f(x0)],'r',x0,stx0,'blackx');
                
                
                axis tight;
                hold off;
                pause;
                if ~ishandle(h)
                    return;
                end
            end
            close(h);
        end
    end
    
end