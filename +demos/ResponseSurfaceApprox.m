classdef ResponseSurfaceApprox < handle
    % This demo shows how the VKOGA algorithm approximates a 2D-1D response
    % surface for different kernel configurations.
    %
    % Just invoke demos.ResponseSurfaceApprox.run
    
properties
    EC;
    XYRange = -5:.1:5;
    PM;
    TrainPerc = .2;
    NumSteps = 10;
    LoopCallback;
end

properties(SetAccess=private)
    X;
    Y;
    Z;
    atd;
end
    
methods
    function this = ResponseSurfaceApprox
        ec = kernels.config.ExpansionConfig;
        ec.Prototype.Kernel = kernels.Wendland;
        ec.StateConfig = kernels.config.WendlandConfig('G',10,'S',2,'Dim',2);
        this.EC = ec;
        pm = PlotManager(false,3,2);
        pm.LeaveOpen = true;
        this.PM = pm;            
    end
    
    function Z = getSurface(this, X, Y)
        Z = sin(X/2).*Y.^2+5*X+Y;
        this.Z = Z;
    end
    
    function alg = setupVKOGA(this)
        alg = approx.algorithms.VKOGA;
        alg.MaxExpansionSize = this.NumSteps;
        alg.UsefPGreedy = 0;
        alg.ExpConfig = this.EC;
    end
    
    function kexp = compute(this)
        this.init;
        alg = this.setupVKOGA;
        kexp = alg.computeApproximation(this.atd);
    end
    
    function init(this, sel)
        [X,Y] = meshgrid(this.XYRange);
        if nargin < 2
            tot = numel(X);
            sel = randperm(tot);
            n = round(this.TrainPerc*tot);
        end
        Z = this.getSurface(X,Y);
        this.Z = Z;
        fxi = Z(sel(1:n));
        xi = [X(sel(1:n)); Y(sel(1:n))];
        atd = data.ApproxTrainData(xi,[],[]);
        atd.fxi = fxi;
        this.X = X;
        this.Y = Y;
        this.atd = atd;
    end
    
    function plotOriginal(this)
        h = this.PM.nextPlot('fun','Original response surface','x_1','x_2');
        surfl(h, this.X, this.Y, this.Z);
        zlabel('f(x)');
        shading interp;
    end
    
    function plotOriginalWithTD(this)
        h = this.PM.nextPlot('fun_atd','Original response surface with training data','x','f(x)');
        surf(h, this.X, this.Y, this.Z);
        shading interp;
        hold(h,'on');
        xi = this.atd.xi.toMemoryMatrix;
        plot3(h,xi(1,:),xi(2,:),fxi,'r.');
    end
    
    function plotTrainData(this)
        h = this.PM.nextPlot('traindata','Training data','x_1','x_2');
        xi = this.atd.xi.toMemoryMatrix;
        plot3(h,xi(1,:),xi(2,:),fxi,'r.');
    end
    
    function iterationPlots(this, kexp, pm)
        if nargin < 3
            pm = this.PM;
        end
        allxi = [this.X(:) this.Y(:)]';

        %         aZ = kexp.evaluate(allxi);
        %         aZ = reshape(aZ,size(X,1),[]);
        %         h = pm.nextPlot('approx','Approx','x','f(x)');
        %         surf(h, X, Y, aZ);
        %         h = pm.nextPlot('error','Error','x','f(x)');
        %         surf(h, X, Y, abs(Z-aZ));
        args = {'EdgeColor','interp','FaceAlpha',.2,'EdgeColor','none','FaceColor','red'};

        centerpos = Utils.findVecInMatrix(allxi,kexp.Centers.xi);

        %% Iteration plots
        zlim = [];
        for s = 0:this.NumSteps
            if s == 0
                fxi = zeros(size(this.X));
            else
                k = kexp.getSubExpansion(s);
                fxiflat = k.evaluate(allxi);
                fxi = reshape(fxiflat,size(this.X,1),[]);
            end
            if isa(pm,'PlotManager')
                h = pm.nextPlot(sprintf('iter%d',s),...
                    sprintf('Approximation at step %d',s),'x_1','x_2');
            else
                h = pm;
                cla(h);
            end
            surf(h, this.X, this.Y, fxi,'EdgeColor','interp');
            zlabel('f(x)');
            hold(h,'on');
            surf(h, this.X, this.Y, this.Z,args{:});
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
            %axis(h,'tight');
            if ~isempty(this.LoopCallback)
                this.LoopCallback(this,s,h);
            end
            if isempty(zlim)
                zlim = get(h,'ZLim');
            else
                set(h,'ZLim',zlim);
            end
        end

        if nargin < 2
            pm.done;
        end
    end
    
end

methods(Static)
    function run
        % This demo shows how the VKOGA algorithm approximates a 2D-1D response
        % surface for different kernel configurations.

        rsa = demos.ResponseSurfaceApprox;
        kexp = rsa.compute;
        rsa.iterationPlots(kexp);
    end
end
end