function GMA_xFeatDemo
%GMA_XFEATDEMO Produces the graphics used in the GMA Salzburg Talk
% @author Daniel Wirtz @date 26.08.2010

Gamma = 2;

bell = kernels.GaussKernel(Gamma);
bell.NewtonTolerance = 1e-3;
h = figure(1);

maxX = 4;
X = linspace(0,maxX,35);

f = @bell.evaluateScalar;
df = @bell.evaluateD1;
ddf = @bell.evaluateD2;

xstart = repmat(bell.x0/2,1,size(X,2));
xfeats = ModifiedNewton(xstart,X);

for k=1:length(X)
    if mod(k,2) == 0
        y = X(k);
        % Get precomputed values
        dxf = df(xfeats);
        
        % Get optimization function
        opt = optFun(X,repmat(y,1,size(X,2)),bell.x0,bell.PenaltyFactor);
        
        % Bound plot
        hlp = (maxX-y)*dxf(k); % also used in plot of tangent!
        opt(opt < hlp) = hlp;
        hlp2 = max(1.5,max(df(X)));
        opt(opt > hlp2) = hlp2;
        
        % Plot!
        ph = plot(X,f(X),'r',X,opt,'b',X,abs(df(X)),'m',X,abs(dxf),'g',[0 y],f(y) + [-y*dxf(k) 0],...
           'black--',y,f(y),'bo',[xfeats(k) xfeats(k) y],[f(xfeats(k)) abs(dxf(k)) abs(dxf(k))],'blackx',X,0,'black');
        legend('Gaussian','Optim. fcn','Absolute Gaussian derivative','Max local derivatives','Max local gradient around y','y','x_y (correlation)');
        set(ph,'LineWidth',2);
        set(ph(end-2),'MarkerSize',5);
        pause;
    end
end

close(h);

    function xtmp = ModifiedNewton(xstart, y)
        xtmp = xstart+2*bell.NewtonTolerance;
        x = xstart;
        
        while any(abs(xtmp-x) > bell.NewtonTolerance)
            
            [g,dg] = optFun(x,y,bell.x0,bell.PenaltyFactor);
            
            xtmp = x;
            x = x - g./dg;
        end
        
    end

    function [g,dg] = optFun(x,y,x0,p)
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

end

