classdef LogNorm
    
    methods(Static)
        
        %% Comparison of similarity transformed jacobian log norms to full log norms
        function [aln, times, st_sizes] = compareSimTransJac_FullJac(m, st_sizes)
            % Computes logarithmic norms of similarity transformed
            % jacobians using the model's offline data, containing `N`
            % trajectory samples.
            %
            % Parameters:
            % m: A BaseFullModel whos offlineGenerations have been run and
            % that uses a DEIMEstimator. @type models.BaseFullModel
            % st_sizes: The sizes `s_1,\ldots,s_n` to use for the
            % similarity transformation. If left empty, all from one to the
            % models ErrorEstimator.JacSimTransMaxSize are used
            % (extensive!) @type rowvec<integer>
            %
            % Return values:
            % aln: A `n \times N` matrix with the approximated logarithmic
            % norms in rows for each sim. trans. size.
            % times: A `n \times N` matrix with the computation time for
            % the logarithmic norm
            % st_sizes: The effectively used sizes (If st_sizes is given as
            % parameter, it is looped through)
            
            % Get required data
            jtd = m.Data.JacobianTrainData;
            QFull = m.ErrorEstimator.QFull;
            f = m.System.f;
            
            % Input checks
            if nargin < 2
                st_sizes = 1:size(QFull,2);
            end
            
            % Preps
            stn = length(st_sizes);
            n = size(jtd.xi,2);
            aln = zeros(stn,n);
            times = zeros(stn,n);
            pi = tools.ProcessIndicator('Computing approximated log norms for %d full jacobians and %d sim. trans. orders',...
                n*stn,false,n,stn);
            for nr = 1:n
                J = f.getStateJacobian(jtd.xi(:,nr),jtd.ti(nr),jtd.mui(:,nr));
                for snr = 1:stn
                    Qk = QFull(:,1:st_sizes(snr));
                    t = tic;
                    aln(snr,nr) = general.Utils.logNorm(Qk'*J*Qk);
                    times(snr,nr) = toc(t);
                    pi.step;
                end
            end
            pi.stop;
            if nargout == 0
                testing.DEIM.compareSimTransJac_FullJac_plots(m, aln, times, st_sizes);
            end
        end
        
        function pm = compareSimTransJac_FullJac_plots(m, aln, times, st_sizes, pm)
            if nargin < 5
                pm = tools.PlotManager(false,2,2);
            end
            jstd = m.Data.JacSimTransData;
            pm.nextPlot('correct_ln',sprintf('Full logarithmic norms\nJacobian full dimension %dx%d',...
                m.System.f.fDim,m.System.f.xDim));
            plot(jstd.LogNorms);
            
            if size(aln,1) > 1
                [X,Y] = meshgrid(1:size(aln,2),st_sizes);
                LN = repmat(jstd.LogNorms,size(aln,1),1);
            end
            
            if size(aln,1) == 1
                h = pm.nextPlot('approx_ln',...
                    sprintf('Approx logarithmic norms\nSimilarity transformation size %d',...
                    st_sizes));
                plot(h, aln);
            else
                h = pm.nextPlot('approx_ln','Approx logarithmic norms',...
                    'training sample','similarity transform size');
                tools.LogPlot.nicesurf(h, X, Y, aln);
            end
            
            h = pm.nextPlot('rel_err','Relative error','training sample',...
                'similarity transform size');
            tools.LogPlot.logsurf(h,X,Y,abs((LN-aln)./LN));
            
            sv = m.ErrorEstimator.QSingVals;
            s = m.ErrorEstimator.JacSimTransMaxSize;
            h = pm.nextPlot('sing_vals',...
                sprintf('Singular values of SVD of eigenvector matrix\nBlack line: Max size of sim. trans., singular value= %g',...
                sv(s)),'POD Modes','Singular values');
            semilogy(h,sv);
            hold on;
            plot(h,[s s+eps],[min(sv) max(sv)],'k');
            hold off;
            
            h = pm.nextPlot('comp_times',...
                sprintf('Average computation times (over %d values)\nfor log norms of similarity transformed matrices',...
                    size(times,2)),'similarity transformation size','mean computation time');
            me = mean(times,2);
            semilogy(h,st_sizes,me,'-s');
            if nargout < 1
                pm.done;
            end
        end
        
        %% Comparison of similarity transformed DEIM-approximated jacobian log norms to full log norms
        function [aln, times, jtimes, deim_orders, st_sizes] = ...
                compareSimTransDEIMJac_FullJac(m, deim_orders, st_sizes)
            % Computes logarithmic norms of similarity transformed AND
            % matrix DEIM approximated jacobians using the model's offline
            % data, containing `N` trajectory samples.
            %
            % Parameters:
            % m: A BaseFullModel whos offlineGenerations have been run and
            % that uses a DEIMEstimator. @type models.BaseFullModel
            % deim_orders: The DEIM orders `d_1,\ldots,d_m` to set for the
            % matrix DEIM of the jacobian. @type rowvec<integer>
            % st_sizes: The sizes `s_1,\ldots,s_n` to use for the
            % similarity transformation. If left empty, all from one to the
            % models ErrorEstimator.JacSimTransMaxSize are used
            % (extensive!) @type rowvec<integer>
            %
            % Return values:
            % aln: A `m \times n \times N` matrix with the approximated logarithmic
            % norms in rows for each sim. trans. size.
            % times: A `m \times n \times N` matrix with the computation times for
            % the logarithmic norm
            % jtimes: A `m \times n \times N` matrix with the computation
            % times for the sim. trans. DEIM approximated jacobians
            % deim_orders: The effectively used DEIM orders (If given as
            % parameter, it is looped through)
            % st_sizes: The effectively used sizes (If given as
            % parameter, it is looped through)
            
            % Get required data
            jtd = m.Data.JacobianTrainData;
            %zi = m.Data.W'*jtd.xi;
            e = m.ErrorEstimator;
            jd = e.JacMDEIM;
            
            % Input checks
            if nargin < 3
                st_sizes = [1:e.JacSimTransMaxSize 0];
                if nargin < 2
                    deim_orders = 1:jd.MaxOrder;
                end
            end
            
            % Preps
            stn = length(st_sizes);
            n = size(jtd.xi,2);
            no = length(deim_orders);
            
            aln = zeros(no,stn,n);
            times = aln;
            jtimes = aln;
            pi = tools.ProcessIndicator('Computing approximated log norms over %d DEIM orders and %d sim. trans. sizes on %d training values',...
                numel(aln),false,no,stn,n);
            for onr = 1:no
                jd.Order = deim_orders(onr);
                for snr = 1:stn
                    jd.setSimilarityTransform(e.QFull(:,1:st_sizes(snr)));
                    for nr = 1:n
                        t = tic;
                        J = jd.evaluate(jtd.xi(:,nr),jtd.ti(nr),jtd.mui(:,nr));
                        jtimes(onr,snr,nr) = toc(t);
                        t = tic;
                        aln(onr,snr,nr) = general.Utils.logNorm(J);
                        times(onr,snr,nr) = toc(t);
                        pi.step;
                    end
                end
            end
            pi.stop;
        end
        
        function pm = compareSimTransDEIMJac_FullJac_plots(m, aln, ...
                times, jtimes, deim_orders, st_sizes, pm)
            % See WSH12 tests_burgers for likely better plot routine
            if nargin < 7
                pm = tools.PlotManager(false,2,2);
            end
            
            jstd = m.Data.JacSimTransData;
            n = size(aln,3);
            ln = reshape(jstd.LogNorms,1,1,[]);
            ln = repmat(ln,[length(deim_orders), length(st_sizes), 1]);
            
            abserr = abs(aln-ln);
            relerr = abs(abserr ./ ln);
            relerr(isnan(relerr)) = 0;
            
            % Aggregate errors over trajectory sample data
            %aln = sqrt(sum(aln.^2,3));
            abserr = mean(abserr,3);
            relerr = mean(relerr,3);
            
            % Mean the times too
            ttimes = times + jtimes;
            times = mean(times,3);
            jtimes = mean(jtimes,3);
            ttimes = mean(ttimes,3);
            
            [X,Y] = meshgrid(st_sizes,deim_orders);
            
            h = pm.nextPlot('abs_err',sprintf('Mean absolute approximation error over %d samples',n),...
                'Similarity transformation size','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,abserr);
            
            h = pm.nextPlot('rel_err',sprintf('Mean relative approximation error over %d samples',n),...
                'Similarity transformation size','DEIM order');
            tools.LogPlot.logsurf(h,X,Y,relerr);
            
            h = pm.nextPlot('comp_times_j',...
                sprintf('Average computation times over %d values\nfor matrix DEIM jacobian evaluation',...
                    n),'Similarity transformation size','DEIM order');
            tools.LogPlot.nicesurf(h,X,Y,jtimes);
            
            h = pm.nextPlot('comp_times',...
                sprintf('Average computation times over %d values\nfor log norm computation of sim.trans. matrix DEIM jacobian',...
                    n),'Similarity transformation size','DEIM order');
            tools.LogPlot.nicesurf(h,X,Y,times);
            
            h = pm.nextPlot('comp_times_total',...
                sprintf('Average total computation times over %d values\n(jac comp + log norm comp)',...
                    n),'Similarity transformation size','DEIM order');
            tools.LogPlot.nicesurf(h,X,Y,ttimes);
            
            if nargout < 1
                pm.done;
            end
        end
        
        function res = getApproxLogNormsAtPos(mo, x, t, mui)
            % Computes logarithmic norms of similarity transformed AND
            % matrix DEIM approximated jacobians at the position given by (x,t) over all
            % given parameters mui.
            %
            % Parameters:
            % m: A BaseFullModel whos offlineGenerations have been run and
            % that uses a DEIMEstimator. @type models.BaseFullModel
            % x: The state space location to use. @type colvec<double>
            % t: The time t that belongs to the state space vector x. @type double
            % mui: A matrix of `\mu` values to compute the approximated logarithmic norms for.
            % If left empty, the model's ParamSamples are used.
            % @type matrix<double>
            %
            % Return values:
            % res: A struct with multiple fields. @type struct
            %   aln: The approximate log norms @type rowvec<double>
            %   times: The computation times for the log norm @type rowvec<double>
            %   jtimes: The computation times for the jacobian @type rowvec<double>
            %   mJ: The MDEIM order used @type integer
            %   k: The similarity transformation size used @type integer
            
            % Get required data
            e = mo.ErrorEstimator;
            if isa(mo,'models.ReducedModel')
                fm = mo.FullModel;
            else
                fm = mo;
            end
            jd = e.JacMDEIM;
            
            % Input checks
            if nargin < 4
                mui = fm.Data.ParamSamples;
            end
            
            % Preps
            n = size(x,2);
            if ~isempty(t)
                if length(t) == 1
                    t = ones(1,n)*t;
                elseif length(t) ~= n
                    error('Size of t vector must match column count of x');
                end
            end
            m = size(mui,2);
            aln = zeros(n,m);
            times = aln;
            jtimes = aln;
            pi = tools.ProcessIndicator('Computing approximated log norms for %d parameter values at %d locations',n*m,false,m,n);
            for i = 1:n
                for j = 1:m
                    tc = tic;
                    J = jd.evaluate(x(:,i),t(i),mui(:,j));
                    jtimes(i,j) = toc(tc);
                    tc = tic;
                    aln(i,j) = general.Utils.logNorm(J);
                    times(i,j) = toc(tc);
                    pi.step;
                end
            end
            pi.stop;
            res.aln = aln;
            res.times = times;
            res.jtimes = jtimes;
            res.mJ = e.JacMatDEIMOrder;
            res.k = e.JacSimTransSize;
            res.mui = mui;
            res.munames = {mo.System.Params(:).Name};
        end
        
        function getApproxLogNormsAtPos_plots(res, pm)
            str = sprintf('\nm_J = %d, k=%d, %d locations',res.mJ, res.k, size(res.aln,1));
            if size(res.mui,1) == 1
                if size(res.aln,1) == 1
                    h = pm.nextPlot('aln_min',['Log norms at given location' str],res.munames{1});
                    tools.LogPlot.cleverPlot(h,res.mui,res.aln,'b');
                else
                    [X,Y] = meshgrid(res.mui,1:size(res.aln,1));
                    h = pm.nextPlot('aln_min',['Log norms over all locations' str],res.munames{1},'location nr');
                    surf(h,X,Y,res.aln,'EdgeColor','none','FaceColor','interp');
                end
            elseif size(res.mui,1) == 2
                tri = delaunay(res.mui(1,:),res.mui(2,:));
                if size(res.aln,1) == 1
                    h = pm.nextPlot('aln_min',['Log norms at given location' str],res.munames{:});
                    trisurf(tri,res.mui(1,:),res.mui(2,:),res.aln,...
                        'EdgeColor','none','FaceColor','interp','Parent',h);
                else
                    %% Animation
%                     h = pm.nextPlot('aln','Approximated Jacobian log norms over sample parameters');
%                     for k=1:size(res.aln,1)
%                         trisurf(tri,res.mui(1,:),res.mui(2,:),res.aln(k,:),...
%                             'EdgeColor','none','FaceColor','interp','Parent',h);
%                         title(h,sprintf('Approximated Jacobian log norms over sample parameters, pos: %d',k));
%                         view(90,0);
%                         pause;
%                     end
                    %% Min/Max plot
                    alnmin = min(res.aln,[],1);
                    alnmax = max(res.aln,[],1);
                    h = pm.nextPlot('aln_min',['Minimum log norms over all locations' str],res.munames{:});
                    trisurf(tri,res.mui(1,:),res.mui(2,:),alnmin,...
                        'EdgeColor','none','FaceColor','interp','Parent',h);
                    h = pm.nextPlot('aln_min',['Max log norms over all locations' str],res.munames{:});
                    trisurf(tri,res.mui(1,:),res.mui(2,:),alnmax,...
                        'EdgeColor','none','FaceColor','interp','Parent',h);
                    h = pm.nextPlot('aln_min',['Max-Min log norms over all locations' str],res.munames{:});
                    trisurf(tri,res.mui(1,:),res.mui(2,:),alnmax-alnmin,...
                        'EdgeColor','none','FaceColor','interp','Parent',h);
                end
            else
                error('Cannot plot results for more than two parameters');
            end
        end
        
        function [res, mScale, MScale, pos, l, sel, seli] = CompLogNorms(m, numt)
            % LogNorm: 
            %
            %
            %
            % @author Daniel Wirtz @date 2012-05-08
            %
            % @new{0,6,dw,2012-05-08} Added this function.
            %
            % This class is part of the framework
            % KerMor - Model Order Reduction using Kernels:
            % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
            % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
            % - \c License @ref licensing

            atd = m.Data.ApproxTrainData;
            N = size(atd.xi,2);

            if nargin < 2
                numt = N;
            else
                numt = min(numt, N);
            end
            seed = 1;
            % Show plots for each numt value
            doplot = 4==5;
            % Sort training points after norm values (ascending)
            dosort = 4==5;
            % Plot negative local logarithmic norms
            plotneg = 4==4;

            r = RandStream('mt19937ar','Seed',seed);

            %sel = unique(r.randi(N,1,num));
            sel = 1:N;
            seli = unique(r.randi(N,1,numt));
            n = length(seli);

            pi = tools.ProcessIndicator('Computing %d local log norms, maxed over %d training points',n,false,n,length(sel));
            ln = zeros(1,n);
            l = ln;
            pos = l;

            if doplot
                pm = tools.PlotManager(false, 2, 2);
            end

            xi = atd.xi(:,sel);
            fxi = atd.fxi(:,sel);
            V = m.Data.V;
            W = m.Data.W;

            %% Compute constant terms
            xinsq = sum(xi.^2);

            % Optional sorting of values
            if doplot && dosort
                [xinsq sortidx] = sort(xinsq);
                xi = xi(:,sortidx);
                fxi = fxi(:,sortidx);
            end

            xifxi = sum(xi.*fxi);
            Vfxi = V'*fxi;
            Vxi = V'*xi;
            VV = 1;%V'*V;
            % if norm(VV) - 1 < 1e-12
            %     VV = 1;
            % end

            % Locally selected "Vz" terms
            zi = W'*atd.xi(:,seli);
            MU = atd.mui(:,seli);
            fVzi = m.System.f.evaluate(V*zi,atd.ti(seli),MU);

            % Only needed for efficient local lipschitz constant computation
            % fxinsq = sum(fxi.*fxi);

            %% Main loop
            for i = 1:n
                z = zi(:,i);
                fz = fVzi(:,i);
                nom = xifxi - z'*Vfxi + z'*(V'*fz);
                nom = nom - fz'*xi;
                denom = xinsq - 2*z'*Vxi + z'*(VV*z);
                Ln = nom ./ denom;
                %denom(zer) = sum((xi(:,zer)-V*z(:,ones(1, length(zer)))).^2);

            %     L = sqrt(abs((fxinsq - 2*fz'*fxi + fz'*fz) ./ denom));
            %     L2 = Norm.L2(fxi - fVzi(:,ones(1,length(sel))*i)) ./ sqrt(denom);

                [~, idx] = sort(denom);
                Ln = Ln(idx);

                zer = denom == 0;
                Ln(zer) = [];

                [ln(i), pos(i)] = max(Ln);

            %     [~, sortidx] = sort(denom);
            %     Ln = Ln(sortidx);
                if doplot
                    ti = sprintf('Idx %d: loc Lip const %g, loc log norm %g',seli(i),l(i),ln(i));
                    ax = pm.nextPlot('blah',ti,'training point','local [lip const (b), log norm (g)]');
                    neg = Ln <= 0;
                    if plotneg
                        Lnn = -Ln;
                    end
                    Ln(neg) = 0;
                    hlp = [L; Ln];
                    if plotneg
                        Lnn(~neg) = 0;
                        hlp = [hlp; Lnn];%#ok
                    end
                    semilogy(ax, hlp');
                end
                pi.step;
            end
            if doplot
                pm.done;
            end
            pi.stop;

            res = data.ApproxTrainData(zi,[],[]);
            res.fxi = ln;        
            res = data.ApproxTrainData.makeUniqueXi(res);
            [res, mScale, MScale] = data.ApproxTrainData.scaleXiZeroOne(res);
        end
    end
end