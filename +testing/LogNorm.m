classdef LogNorm
    
    methods(Static)
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