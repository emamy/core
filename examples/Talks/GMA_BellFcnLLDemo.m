function GMA_BellFcnLLDemo
% Creates the graphics used in the Salzburg GMA Talk.
% @author Daniel Wirtz @date 26.08.2010

h = figure(1);
dt = 0.05;
b = 6;

x0 = sqrt(b/2)+2;
Cint = [2.6 Inf];

x = 0:dt:8.5;
x = union(x,x0);

K = kernels.GaussKernel(b);

f = @(x)K.evaluateScalar(x);
df = @(x)K.evaluateD1(x);
ddf = @(x)K.evaluateD2(x);

fx = f(x);
maxR = K.x0;
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


