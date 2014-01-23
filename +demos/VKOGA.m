classdef VKOGA
% VKOGA: Contains some demo functions for the VKOGA algorithm.
%
% @author Daniel Wirtz @date 2013-01-15
%
% @new{0,7,dw,2013-01-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
        
    methods(Static)
        
        function res = VKOGA_1D_nD(n,fPGreedy,nG)
            % Starts a demo of the approx.algorithms.VKOGA algorithm
            %
            % The functions all live on `[-5,5]` for simplicity and can be
            % of output dimension `1-4`.
            %
            % Parameters:
            % n: The output space dimension. Between 1 and 4. @type integer @default 1
            % fPGreedy: Flag to use f/P-Greedy instead of f-Greedy. @type
            % logical @default false
            % nG: The number of different `\gamma` values to use for the
            % Gaussian. @type integer @default 3
            %
            % Return values:
            % res: A struct with the fields {kexp, atd, alg} for further
            % processing. @type struct
            
            %% Data setup
            x = -5:.1:5;
            fx = [sin(x)+.2*x; cos(x)-.1*x; exp(-x.^2).*sin(x); cos(2*x).*sin(x)];
            if nargin < 3
                nG = 3;
                if nargin < 2
                    fPGreedy = false;
                    if nargin < 1
                         n = 1;
                    elseif n > 2
                        n = size(fx,1);
                    end
                end
            end
            fx = fx(1:n,:);
            atd = data.ApproxTrainData(x,[],[]);
            atd.fxi = fx;
            
            %% Algorithm setup
            alg = approx.algorithms.VKOGA;
            alg.MaxExpansionSize = 50;
            alg.MaxRelErr = 1e-5;
            alg.UsefPGreedy = fPGreedy;
            ec = kernels.config.ExpansionConfig;            
            ec.Prototype.Kernel = kernels.GaussKernel(.8);
            ec.StateConfig = kernels.config.GaussConfig('G',linspace(.5,2,nG));
            alg.ExpConfig = ec;

            kexp = alg.computeApproximation(atd);
            
            %% Plot approximated function
            m = length(alg.Used);
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            
            h1 = pm.nextPlot('fun','Function','x','f(x)');
            plot(h1,x,[fx; kexp.evaluate(x)]');
            
            h2 = pm.nextPlot('nfun','Newton Basis Functions on training data','x','N(x)');
            plot(h2,x,alg.bestNewtonBasisValuesOnATD); 
            
            h3 = pm.nextPlot('err','Absolute error','x','|f(x)-f^m(x)|');
            ph = LogPlot.cleverPlot(h3,1:m,alg.MaxErrors(:,1:m));
            set(ph(alg.ExpConfig.vBestConfigIndex),'LineWidth',2);
            
            h4 = pm.nextPlot('relerr','Relative error','x','|(f(x)-f^m(x))/f(x)|');
            ph = LogPlot.cleverPlot(h4,1:m,alg.MaxRelErrors(:,1:m));
            set(ph(alg.ExpConfig.vBestConfigIndex),'LineWidth',2);
            pm.done;
            
            res.kexp = kexp;
            res.atd = atd;
            res.alg = alg;
        end
        
        function IterationPlots(res, steps, pm)
            % Demonstrates the VKOGA iterations during approximation
            % computations.
            %
            % Parameters:
            % res: The output of demos.VKOGA#demo_VKOGA_1D_nD @type struct
            % steps: How many steps to illustrate @type integer @default 4
            % pm: A PlotManager instance to use for plots. @type
            % PlotManager @default 3x3 Subfigures
            if nargin < 3
                pm = PlotManager(false,3,3);
                pm.LeaveOpen = true;
                if nargin < 2
                    steps = 9;
                end
            end
            
            if res.alg.UsefPGreedy
                fprintf(2,'Warning: f/P-Greedy selection criteria was used. Error plots are for f-Greedy case.\n');
            end
            ms = 7;
            fx = res.atd.fxi.toMemoryMatrix;
            x = res.atd.xi.toMemoryMatrix;
            
            %% Zero function plot
            k = [];
            h0 = doPlot(0,zeros(size(fx)));
            
            %% Iteration plots
            for s = 1:steps
                k = res.kexp.getSubExpansion(s);
                fxi = k.evaluate(x);
                hn = doPlot(s,fxi);
                axis(hn,axis(h0));
            end
            
            if nargin < 3
                pm.done;
            end
            
            function h1 = doPlot(s,fxi)
                h1 = pm.nextPlot(sprintf('step%d',s),...
                    sprintf('Iteration %d',s),'x','f(x)');
                plot(h1,x,fx','LineWidth',2); 
                hold(h1,'on');
                
                err = fx-fxi;
                allerr = Norm.L2(err);
                plot(h1,x,allerr,'r'); 
                plot(h1,x,fxi','--'); 
                if ~isempty(k)
                    plot(h1,k.Centers.xi,k.evaluate(k.Centers.xi)','k.','MarkerSize',17);
                end
                
                % Plot (next) max errors
                [v, idx] = max(allerr);
                plot(h1,x(idx),v,'ro','MarkerSize',ms); 
                [~, idx] = max(abs(err),[],2);
                if s == 2
                    for l=1:length(idx)
                        plot(h1,[x(idx(l)) x(idx(l))+eps],[fx(l,idx(l)) fxi(l,idx(l))],'k--');
                        plot(h1,[x(idx(l)) x(idx(l))+eps],[fx(l,idx(l)) fxi(l,idx(l))],'kx','MarkerSize',ms);
                    end
                end
                axis(h1,'tight');
                axis(h1,axis(h1)*1.03);
            end
        end
        
        function NewtonBasis_Schaback
            % The demo of the schaback paper @cite PS11 for the
            % function-dependent Newton basis.
            %
            % Original source code from the link below, adopted to fit into
            % KerMor.
            %
            % See also:
            % \li http://num.math.uni-goettingen.de/schaback/research/papers/BfKBS.tgz
            % \li http://num.math.uni-goettingen.de/schaback/research/group.html
            
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            
            % get cut circle
            [X, ind, Xemesh, Yemesh] = cutcircle(.05);
            
            % assemble training data
            atd = data.ApproxTrainData(X,[],[]);
            fxi = exp(abs(X(1,:)-X(2,:)))-1;
            atd.fxi = fxi;
            
            %% Algorithm setup
            kexp = kernels.KernelExpansion;
            
            alg = approx.algorithms.VKOGA;
            % Maximal number of data sites to be finally used
            alg.MaxExpansionSize = 200;
            % tolerance for residual of function recovery
            alg.MaxAbsResidualErr = 1e-6;
            % "disable" maxrelerr check by ridiculous small amount
            alg.MaxRelErr = 1e-10;
            alg.UsefPGreedy = false;
            ec = kernels.config.ExpansionConfig;
            
            % Gaussian
%             kexp.Kernel = kernels.GaussKernel;
%             ec.StateConfig = kernels.config.RBFConfig('G',[1 2 3]); 
            
            % Wendland functions
            kexp.Kernel = kernels.Wendland;
            wc = kernels.config.WendlandConfig('G',[1 2 3],'S',[1 1 1]);
            wc.Dimension = 2;
            ec.StateConfig = wc;
            
            alg.ExpConfig = ec;

            alg.computeApproximation(kexp, atd);
            
            h = pm.nextPlot('circ','Selected data points');
            plot(h,X(1,:),X(2,:),'.',X(1,alg.Used),X(2,alg.Used),'ro');
            
            Z = zeros(size(Xemesh));
            Z(ind) = fxi;
            h = pm.nextPlot('fun','Original function');
            surf(h,Xemesh,Yemesh,Z);
            
            h = pm.nextPlot('afun','Approximation');
            aZ = Z;
            aZ(ind) = kexp.evaluate(X);
            surf(h,Xemesh,Yemesh,aZ);
            
            h = pm.nextPlot('err','Error');
            surf(h,Xemesh,Yemesh,abs(Z-aZ));
            
            for n=1:8
                plotBasis(n);
            end
            pm.done;
            
            function plotBasis(idx)
                h = pm.nextPlot(sprintf('bfun%d',idx),sprintf('Newton basis function %d',idx),'x','y');
                nv = zeros(size(Xemesh));
                nv(ind) = alg.bestNewtonBasisValuesOnATD(:,idx);
                surf(h,Xemesh,Yemesh,nv);
            end
            
            function [Xe, ind, xm, ym]=cutcircle(he)
                % generates regular points Xe in circle with lower left quadrangle cut out
                % The index list ind indicates the positions in meshgrid (xm,ym) 
                % when the meshgrid is written in a single point list (xee, yee)
                dx=2; % space dimension
                dy=dx;
                % he=0.01; % 1 D stepsize for evaluation grid. 
                % The grid can consist of very many points.
                % he=0.01 gives N=40401 points to work on
                [xm, ym]=meshgrid(-1:he:1,-1:he:1);
                xee=xm(:);
                yee=ym(:);
                ind=find((xee.^2+yee.^2<=1)&((xee>=0)|(yee>=0)) );
%                 msk=zeros(size(xee));
%                 msk(ind)=1.0;
                Xe=[xee(ind), yee(ind)]';
            end
        end
    end
    
    
    
end