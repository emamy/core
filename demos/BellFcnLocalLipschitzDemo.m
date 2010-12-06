function BellFcnLocalLipschitzDemo( x0, C )
%ERRORESTDEMO Demo for the monotone radial basis functions error estimator.

h = figure(1);
dt = 0.05;
b = 6;

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

K = kernels.GaussKernel(b);

f = @(x)K.evaluateScalar(x);
df = @(x)K.evaluateD1(x);
ddf = @(x)K.evaluateD2(x);

fx = f(x);
maxR = K.x0;

% precompute xfeats-vector
xfeats = newton(maxR+sign(maxR-x),x,f,df,ddf,...
    K.NewtonTolerance,K.x0,K.PenaltyFactor);

for C = Cint
    maxder = ones(size(x))*abs(df(K.x0));
    e1 = maxder;
    e2 = maxder;
    e3 = maxder;
    minder = zeros(size(x));
    maxd = zeros(size(x));
    for i = 1:length(x)
        d = abs(x(i));

        % Less efficient, but shows advantage of estimation on [maxR, xfeat]
        if d - C - maxR > 0
            e1(i) = abs(df(d-C));
            e2(i) = (f(d-C) - f(d)) / C;
            e3(i) = e2(i);
        elseif d + C - maxR < 0
            e1(i) = abs(df(d+C));
            e2(i) = (f(d) - f(d+C)) / C;
            e3(i) = e2(i);
        else
            xfeat = newton(maxR+sign(maxR-d),d,f,df,ddf,...
               K.NewtonTolerance,K.x0,K.PenaltyFactor);
            if d + C - xfeat < 0
                e3(i) = (f(d) - f(d+C)) / C;
            elseif d - C - xfeat > 0
                e3(i) = (f(d-C) - f(d)) / C;
            else
                e3(i) = abs(df(xfeat));
            end
        end
        
        
        % More efficient method (even improved e2, e1)
%         xfeat = maxR;
%         if abs(d-maxR) < C
%             xfeat = newton(maxR+sign(maxR-d),d,f,df,ddf,...
%                 K.NewtonTolerance,K.x0,K.PenaltyFactor);
%         end
%         if d + C - xfeat < 0
%             e1(i) = abs(df(d+C));
%             e2(i) = (f(d) - f(d+C)) / C;
%             e3(i) = e2(i);
%         elseif d - C - xfeat > 0
%             e1(i) = abs(df(d-C));
%             e2(i) = (f(d-C) - f(d)) / C;
%             e3(i) = e2(i);
%         else
%             e3(i) = abs(df(xfeat));
%         end
        
        
        % Minimum possible derivative
        di = abs(x - x(i));
        idx = di <= C;
        minder(i) = max(abs(fx(idx) - exp(-x(i).^2/b)) ./ di(idx));
        
        % maxD-Comp
        d = (x(i)-x0);
        maxd(i) = (fx(i)-exp(-x0.^2/b))/d;
        
    end
    
    ph = plot(x,f(x),'r',x,abs(df(x)),'r--',x,e1,'b',x,e2,'g',x,e3,'m');
    set(ph,'LineWidth',2);
    hold on;
    %plot(x,e12,'b--',x,e22,'g--',x,e32,'m--');
    title(sprintf('C=%f',C));
    
    %legend('f(x)','Df(x)','Local Gradient with x_0','Local Secant with x_0','Local secant with x_y');
    
    skip = round(1/(3*dt));
    ph = plot(x(1:skip:end),minder(1:skip:end),'b^');
    set(ph,'MarkerSize',6);
    %plot(x,minder,'r');
    %legend('f(x)','Df(x)','Local Gradient with x_0','Local Secant with x_0','Local secant with x_y','Minimal Lipschitz constant');
    legend('Gaussian','Abs. Gaussian derivative','Local Gradient with x_0','Local Secant with x_0','Local secant with x_y','Minimal Lipschitz constant');
    
    % C-range plot
    %     plot([x0-C x0-C+eps], [0 1], 'black--',[x0 x0+eps], [0 1], 'black',[x0+C x0+C+eps], [0 1], 'black--');
    
    % 0-C plot
    %plot([C C+eps],[0 1],'black.-.');
    
    % penalty-fkt
    %     new = df(x) - (f(x)-f(x0))./(x-x0) - max(0,sign(x0-maxR)*(x-maxR)).^2;
    %     new(abs(new) > 1) = -1;
    %     plot(x,new,'b');
    %plot(x,abs(df(xfeats)),'g--');
    
    % maxR-Stelle (sqrt(b/2))
    %plot(x,0,'black',[maxR maxR+eps],[0 1],'y');
    
    % maxderivative at x0-Plot
    %plot(x,maxd,'black');
    %plot(x,max(maxd),'black');
    
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
end
close(h);

    function xtmp = newton(x, xfix, f, df, ddf, tol, x0, p)
        xtmp = x+2*tol;
        while any(abs(xtmp-x) > tol)
            [g,dg] = kernels.BellFunction.optFun(x,xfix,f,df,ddf,x0,p);
            xtmp = x;
            x = x - g./dg;
        end
    end
end


