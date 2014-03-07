function ResponseSurfaceApprox(ec, pm)
% This demo shows how the VKOGA algorithm approximates a 2D-1D response
% surface for different kernel configurations.
%
% Parameters:
% ec: An ExpansionConfig instance for approximation configuration.
% @type kernels.config.ExpansionConfig
% pm: A PlotManager instance to use @type PlotManager @default
% PlotManager(false,2,2)

if nargin < 2
    pm = PlotManager(false,2,2);
    pm.LeaveOpen = true;
    if nargin < 1
        ec = kernels.config.ExpansionConfig;
        ec.Prototype.Kernel = kernels.Wendland;
        ec.StateConfig = kernels.config.WendlandConfig('G',10,'S',2,'Dim',2);
    end
end

%% Data setup
[X,Y] = meshgrid(-5:.1:5);
Z = sin(X/2).*Y.^2+5*X+Y;
tot = numel(X);
sel = randperm(tot);

n = round(.2*tot);
fxi = Z(sel(1:n));
xi = [X(sel(1:n)); Y(sel(1:n))];
allxi = [X(:) Y(:)]';

h = pm.nextPlot('fun','Original response surface','x_1','x_2');
surfl(h, X, Y, Z);
zlabel('f(x)');
shading interp;

h = pm.nextPlot('fun_atd','Original response surface with training data','x','f(x)');
surf(h, X, Y, Z);
shading interp;
hold(h,'on');
plot3(h,xi(1,:),xi(2,:),fxi,'r.');

h = pm.nextPlot('traindata','Training data','x_1','x_2');
plot3(h,xi(1,:),xi(2,:),fxi,'r.');
zlabel('f(x)');

atd = data.ApproxTrainData(xi,[],[]);
atd.fxi = fxi;

%% Algorithm setup
alg = approx.algorithms.VKOGA;
alg.MaxExpansionSize = 20;
alg.UsefPGreedy = 0;
alg.ExpConfig = ec;
kexp = alg.computeApproximation(atd);
%         aZ = kexp.evaluate(allxi);
%         aZ = reshape(aZ,size(X,1),[]);
%         h = pm.nextPlot('approx','Approx','x','f(x)');
%         surf(h, X, Y, aZ);
%         h = pm.nextPlot('error','Error','x','f(x)');
%         surf(h, X, Y, abs(Z-aZ));
args = {'EdgeColor','interp','FaceAlpha',.2,'EdgeColor','none','FaceColor','red'};

centerpos = Utils.findVecInMatrix(allxi,kexp.Centers.xi);

%% Iteration plots
for s = 0:10
    if s == 0
        fxi = zeros(size(X));
    else
        k = kexp.getSubExpansion(s);
        fxiflat = k.evaluate(allxi);
        fxi = reshape(fxiflat,size(X,1),[]);
    end
    h = pm.nextPlot(sprintf('iter%d',s),...
        sprintf('Approximation at step %d',s),'x_1','x_2');
    surf(h, X, Y, fxi,'EdgeColor','interp');
    zlabel('f(x)');
    hold(h,'on');
    surf(h, X, Y, Z,args{:});
    if s > 0
        plot3(allxi(1,centerpos(s)),allxi(2,centerpos(s)),...
        fxiflat(1,centerpos(s)),'y*','MarkerSize',20);
        plot3(allxi(1,centerpos(s)),allxi(2,centerpos(s)),...
            fxiflat(1,centerpos(s)),'k.','MarkerSize',20);
        if s > 1
            plot3(allxi(1,centerpos(1:s-1)),allxi(2,centerpos(1:s-1)),...
                fxiflat(1,centerpos(1:s-1)),'r.','MarkerSize',16);
        end
    end
    %     h = pm.nextPlot(sprintf('err_iter%d',s),...
    %         sprintf('Error at step %d',s),'x','f(x)');
    %     surf(h, X, Y, abs(Z-fxi),'EdgeColor','interp');
    axis(h,'tight');
end

if nargin < 2
    pm.done;
end

end