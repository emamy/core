classdef ScalarEpsSVR_SMO < general.regression.BaseScalarSVR
% ScalarEpsSVR_SMO: 
%
%
%
% @author Daniel Wirtz @date 2011-08-09
%
% @new{0,5,dw,2011-08-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        Eps = .1;
        
        StopEps = 1e-3;
                
        Vis = 1;
        
        Version = 1;
    end
    
    methods
        function ai = regress(this, fxi)
            % Make sure it's a row vector.
            fxi = reshape(fxi,1,[]);
            if this.Version == 1
                % Call 1D-Optimizer
                ai = this.regress1D(fxi);
            else
                % Call 2D-Optimizer
                ai = this.regress2D(fxi);
            end
        end
    end
    
    methods(Access=private)
        
        function ai = regress1D(this, fxi)
            
            %% Preps
            this.Vis = 0;
            this.Vis = 1;
%             this.Vis = 2;
            maxcnt = 20000;
            cnt = 1;
            
            n = length(fxi);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            Sall = zeros(1,maxcnt);
            S = n*this.C;
            
            %% Init - cold start
            [a, dW, T] = I0_ColdStart(this, fxi);
            
            %% Init - kernel rule
%             [a, dW, T] = I1_ColdStartKernelRule(this, fxi);
            
            if this.Vis > 0
                if this.Vis > 1
                    h = figure(1);
                end
                err = Sall; 
                Tall = Sall; 
                Eall = Sall;
                DualAll = Sall;
            end
            
            stop = this.StopEps/(2*this.Lambda);
            while S > stop && cnt < maxcnt
                
                if this.Vis > 0
                    afxi = (a(1:n)-a(n+1:end))*this.K;
                    etmp = abs(fxi-afxi);
                    err(cnt) = sum(etmp .* (etmp > this.Eps));
                    if this.Vis > 1
                        figure(h);
                        subplot(1,2,1);
                        plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b');
                        title(sprintf('Current approximation, error: %e',err(cnt)));
                        legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                        axis tight;
                    end

                    if sum(a(1:n).*a(n+1:end)) ~= 0
                        warning('some:id','Invalid ai config. please check now!');
                        keyboard;
                    end
                end
                
                %% Select working set
                [idx, M] = this.WSS0_1D_GetMaxGainIndex(a, dW, n);
                
                %[idx, M] = this.WSS1_1D_GetMaxGradientIndex(a, dW);

                r = max(-a(idx),min(this.C-a(idx),dW(idx)));
                
                %% Misc vars
                % For index: Subtract n if index is for an \alpha^- (linear indexing)
                idx1 = idx - n*(idx > n);
                
                if this.Vis > 1
                    % Compute some stuff again for visual
                    ind = 1:n;
                    ch = [ind(a(n+1:end) == 0) ind(a(1:n) == 0)+n];
                    
                    subplot(1,2,2);
                    r2 = zeros(1,2*n);
                    r2(ch) = max(-a(ch),min(this.C-a(ch),dW(ch)));
                    plot(r2,'b');
                    hold on;
                    g2 = zeros(1,2*n);
                    g2(ch) = r2(ch) .* (dW(ch) - .5*r2(ch));
                    plot(g2,'g');
                    plot(idx,g2(idx),'r.','MarkerSize',6);
                    plot([n n+eps],[0, max([r2 g2])],'black');
                    hold off;
                    title(sprintf('Best gain index: %d, gain: %e, \\alpha_{%d} change: %e',idx, M, idx, r2(idx)));
                    legend('\alpha difference','gain');
                    axis tight;
                    subplot(1,2,1);
                    hold on;
                    plot(idx1,afxi(idx1),'.','MarkerSize',5);
                    hold off;
                end
                
                %% update ai at max gain index
                a(idx) = a(idx) + r;                  
                
                %% Get updated stopping criteria
                T = T + r*(r - 2*dW(idx) + sgn(idx)*fxi(idx1) - this.Eps);

                E = this.C * sum(max(0,min(2-this.Eps,dW)));
                
                %% Update dW - only on side with change
                Ki = this.K(idx1, :);
                dW = dW + sign(idx-n-.5) * r * [Ki -Ki];
                
                %% Duality gap P-W (testing)
                dif = (a(1:n)-a(n+1:end)) * this.K;
                hlp = abs(fxi - dif) - this.Eps;
                hlp(hlp < 0) = 0;
                t_tmp = (dif-fxi)*(a(1:n)-a(n+1:end))' + this.Eps*sum(a);
                e_tmp = this.C * sum(hlp);
                DualAll(cnt) = t_tmp + e_tmp;
                
                S = T + E;
                Sall(cnt) = S;
                Tall(cnt) = T;
                Eall(cnt) = E;
                cnt = cnt+1;
            end
            
            fprintf('Finished after %d/%d iterations.\n',cnt,maxcnt);
            if this.Vis > 0
                figure;
                semilogy(1:cnt,Sall(1:cnt),'r',1:cnt,abs(Tall(1:cnt)),'b',1:cnt,Eall(1:cnt),'g',1:cnt,ones(1,cnt)*stop,'black');
                axis tight;
                title('S, T, E and stopping values'); legend('S = T+E','T','E','stop cond');

                figure;
                semilogy(err(1:cnt));
                axis tight;
                title('Approximation error (by eps-insensitive loss fcn)');
                
                figure;
                semilogy(DualAll(1:cnt));
                axis tight;
                title('Duality Gap');
            end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (a(1:n)-a(n+1:end))';
        end
        
        function ai = regress2D(this, fxi)
            
            %% Preps
            this.Vis = 0;
            this.Vis = 1;
%             this.Vis = 2;
            maxcnt = 10000;
            cnt = 1;
            
            n = length(fxi);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            Sall = zeros(1,maxcnt);
            S = n*this.C;

            %% Init - cold start
            [a, dW, T] = I0_ColdStart(this, fxi);
            
            %% Init - kernel rule
%             [a, dW, T] = I1_ColdStartKernelRule(this, fxi);
            
            if this.Vis > 0
                if this.Vis > 1
                    h = figure(1);
                end
                err = Sall; 
                Tall = Sall; 
                Eall = Sall;
                DualAll = Sall;
            end

            % Initialize i to -1 if MaxGainPlusLast is used.
            i = -1;
            
            stop = this.StopEps/(2*this.Lambda);
            while S > stop && cnt < maxcnt
                
                %% Max gain computation - O(n^2) strategy
                %[i, j] = this.WSS0_BestGainIndices(a, dW, sgn);
                
                %[i, j] = this.WSS1024_RandomAdmissibleIndices(a);
                
                if cnt == 17
                    asdg = 4;
                end
                %[i, j] = this.WSS1_1DMaxGainPlusLast(a, dW, n, i);
                
                [i, j] = this.WSS2_TwoSetMaxGain(a, dW, n);
                
                %[i, j] = this.WSS4_TwoMaxGradients(a, dW);
                
                i1 = i - (i > n)*n;
                j1 = j - (j > n)*n;
                
                if this.Vis > 0
                    afxi = (a(1:n)-a(n+1:end))*this.K;
                    etmp = abs(fxi-afxi);
                    err(cnt) = sum(etmp .* (etmp > this.Eps));
                    if this.Vis > 1
                        figure(h);
                        subplot(1,2,1);
                        plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b');
                        hold on;
                        plot([i1 j1],afxi([i1 j1]),'black.','MarkerSize',7);
                        hold off;
                        title(sprintf('Current approximation, error: %e',err(cnt)));
                        %legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                        axis tight;
                    end
                    
                    %fprintf('%f ',sum(ai(1:n).*ai(n+1:end)));
                    if sum(a(1:n) .* a(n+1:end)) ~= 0
                        warning('some:id','Invalid ai config. please check now!');
                        keyboard;
                    end
                end
                
                %% Compute updates
                kij = this.K(i1,j1);
                % Compute plusminus sign of K_ij update.
                % If difference is smaller than n, they are ++ or -- updates, meaning subtraction.
                % otherwise, is |i-j| \geq n, we have +- or -+ cases, meaning addition.
                % subtraction of .5 avoids having sign(0) = 0.
                pm = sign(abs(i-j) - n - .5);
                r = (dW(i) + pm * kij * dW(j)) / (1-kij^2);
                s = (dW(j) + pm * kij * dW(i)) / (1-kij^2);
                
                % Handle constraint cases
                r2 = []; s2 = [];
                if r > this.C-a(i) || r < -a(i)
                    r = min(this.C-a(i), max(-a(i), r));
                    s2 = min(this.C-a(j), max(-a(j), dW(j) - r * kij));
                end
                if s > this.C-a(j) || s < -a(j)
                    s = min(this.C-a(j), max(-a(j), s));
                    r2 = min(this.C-a(i), max(-a(i), dW(i) - s * kij));
                end
                % Both updates are constrained, search for max gain
                if ~isempty(r2) && ~isempty(s2)
                    % Gain if r is constrained and s accordingly updated
                    gr = r*(dW(i)-.5*r) + s2*(dW(j)-.5*s2) + pm * r * s2 * kij;
                    % Gain if s is constrained and r accordingly updated
                    gs = r2*(dW(i)-.5*r2) + s*(dW(j)-.5*s) + pm * r2 * s * kij;
                    if gr < gs
                        r = r2;
                    else
                        s = s2;
                    end
                end
               
                if this.Vis > 1
                   fprintf('alpha_{%d} change: %e, alpha_{%d} change: %e\n',i,r,j,s);
                end
                
                a(i) = a(i) + r;
                a(j) = a(j) + s;
                
                if any(a < 0) || any(a > this.C)
                    warning('some:id','constraint violation. please check.');
                    keyboard;
                end
                
                %% Update dW - only on side with change
                Ki = this.K(i1, :);
                Kj = this.K(j1, :);
                
                % Update dW
                dW = dW + sign(i-n-.5) * r * [Ki -Ki] ...
                        + sign(j-n-.5) * s * [Kj -Kj];

                %% update stopping crits
                T = T + r*(r - 2*dW(i) - this.Eps + sgn(i)*fxi(i1)) ...
                      + s*(s - 2*dW(j) - this.Eps + sgn(j)*fxi(j1)) ...
                      + pm*r*s*kij;
                  
                %T = T + r*(r - 2*dW(idx) + sgn(idx)*fxi(idx1) - this.Eps);

                E = this.C * sum(max(0,min(2-this.Eps,dW)));
                
                %% Duality gap P-W (testing)
                dif = (a(1:n)-a(n+1:end)) * this.K;
                hlp = abs(fxi - dif) - this.Eps;
                hlp(hlp < 0) = 0;
                t_tmp = (dif-fxi)*(a(1:n)-a(n+1:end))' + this.Eps*sum(a);
                e_tmp = this.C * sum(hlp);
                DualAll(cnt) = t_tmp + e_tmp;

                S = T + E;
                Sall(cnt) = S;
                Tall(cnt) = T;
                Eall(cnt) = E;
                cnt = cnt+1;
            end
            
            fprintf('Finished after %d/%d iterations.\n',cnt,maxcnt);
            if this.Vis > 0
                figure;
                semilogy(1:cnt,Sall(1:cnt),'r',1:cnt,abs(Tall(1:cnt)),'b',1:cnt,Eall(1:cnt),'g',1:cnt,ones(1,cnt)*stop,'black');
                axis tight;
                title('S, T, E and stopping values'); legend('S = T+E','T','E','stop cond');

                figure;
                semilogy(err(1:cnt));
                axis tight;
                title('Approximation error (by eps-insensitive loss fcn)');
                
                figure;
                semilogy(DualAll(1:cnt));
                axis tight;
                title('Duality Gap');
            end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (a(1:n)-a(n+1:end))';
            
        end
        
        function [i, M] = WSS0_1D_GetMaxGainIndex(this, a, dW, n)
            ind = 1:n;
            % Check which alpha values might be changed (only the ones with their partner
            % variable equal to zero)
            ch = [ind(a(n+1:end) == 0) ind(a(1:n) == 0)+n];

            % Find optima for changeable alphas
            r = max(-a(ch),min(this.C-a(ch),dW(ch)));

            % Get gain
            g = r .* (dW(ch) - .5*r);
    
            [M, idxch] = max(g);
            i = ch(idxch);
        end
        
        function [i, M] = WSS1_1D_GetMaxGradientIndex(this, a, dW)
             % S2: Max gradient value for nonbounded \alpha
             [M, i] = max(max([dW .* (a < this.C); -dW .* (a > 0)],[],1));
        end
        
        function [i, j] = WSS1024_RandomAdmissibleIndices(this, a)%#ok
            n = numel(a)/2;
            ind = 1:n;
            ch = [ind(a(n+1:end) == 0) ind(a(1:n) == 0)+n];
            i = ch(randi(length(ch)));
            j = ch(randi(length(ch)));
            while j==i
                j = ch(randi(length(ch)));
            end
        end
        
        function [i, j] = WSS0_BestGainIndices(this, a, dW, sgn)
            % Implements the O(n^2)-WSS0 working set selection strategy.
            
            % Only consider changeable \alpha^{+,-}, i.e. the ones whos partner \alpha equals zero.
            n = numel(a)/2;
            ind = 1:n;
            apch = ind(a(n+1:end) == 0);
            amch = ind(a(1:n) == 0);
            ch = [apch amch];
            Kch = this.K(ch,ch);
            ch(numel(apch)+1:end) = ch(numel(apch)+1:end)+n;
            [wj, wi] = meshgrid(dW(ch));
            ai = repmat(a(ch)', 1, length(ch));
            aj = ai';
            sig = sgn(ch)' * sgn(ch);
            div = 1 ./ (1-Kch.^2);
            div(isinf(div)) = 0;
            
            % Compute update r
            r = (wi - sig .* Kch .* wj) .* div;
            s = r';
            
            %% Handle constraints
            ub = r > this.C - ai;
            lb = r < - ai;
                        
            % r is constrained only
            ronly = (ub | lb) & ~(ub' | lb');
            r(ronly) = min(this.C - ai(ronly),max(- ai(ronly), r(ronly)));
            s(ronly) = min(this.C-aj(ronly), max(-aj(ronly), wj(ronly) - r(ronly) .* Kch(ronly)));
            
            % s is constrained only
            sonly = (ub' | lb') & ~(ub | lb);
            s(sonly) = min(this.C - aj(sonly),max(- aj(sonly), s(sonly)));
            r(sonly) = min(this.C-ai(sonly), max(-ai(sonly), wi(sonly) - s(sonly) .* Kch(sonly)));
            
            % Case of both vars constrained
            % Update s2 and r2 for the cases of the other one is constrained
            both = (ub | lb) & (ub' | lb');
            r1 = zeros(size(r)); r2=r1; s1=r1; s2=r1;
            r1(both) = min(this.C - ai(both),max(- ai(both), r(both)));
            s1(both) = min(this.C-aj(both), max(-aj(both), wj(both) - r1(both) .* Kch(both)));
            s2(both) = min(this.C - aj(both),max(- aj(both), s(both)));
            r2(both) = min(this.C-ai(both), max(-ai(both), wi(both) - s2(both) .* Kch(both)));
            
            % Determine update via get higher gain
            % Gain if r is constrained and s accordingly updated
            g1 = zeros(size(r)); g2 = g1;
            %gr(hlp) = r(hlp) .* (wi(hlp) -.5*r(hlp)) + s2(hlp).*(wj(hlp) -.5 * s2(hlp))- r(hlp) .* s2(hlp).* Kch(hlp);
            g1(both) = r1(both) .* (wi(both) -.5*r1(both)) + s1(both).*(wj(both) -.5 * s1(both))- r1(both) .* s1(both).* Kch(both);
            % Gain if s is constrained and r accordingly updated
            g2(both) = r2(both).* (wi(both) -.5*r2(both))+ s2(both) .*(wj(both) -.5 * s2(both)) - r2(both).* s2(both) .* Kch(both);
            
            % See which gain is higher.
            cmp = false(size(r));
            cmp(both) = g1(both) > g2(both);
            % If gs is larger, r must be updated with the r2 values at those points
            r(both & cmp) = r1(both & cmp);
            s(both & cmp) = s1(both & cmp);
            % If gr is larger, s must be updated with the s2 values at those points
            r(both & ~cmp) = r2(both & ~cmp);
            s(both & ~cmp) = s2(both & ~cmp);
            
            % Compute gain for r_i,s_j combinations
            g = r .* (wi - .5 * r) + s .* (wj - .5 * s) - r .* s .* Kch;
            
            % Extract i,j indices
            [hlp, i] = max(g,[],1);
            [maxg, j] = max(hlp);
            i = i(j);
            
            % Transfer back to full indices
            i = ch(i);
            j = ch(j);
            
            if this.Vis > 1
                [X,Y] = meshgrid(1:2*n,1:2*n);
                subplot(1,2,2);
                G = zeros(2*n,2*n);
                G(ch,ch) = g;
                surf(X,Y,G,'EdgeColor','none');
                hold on;
                plot3(X(i,j),Y(i,j),G(i,j),'black.','MarkerSize',8);
                hold off;
                title(sprintf('Best gain index: (%d,%d), gain: %e',i,j, maxg));
                axis tight;
            end
        end
        
        function [i,j] = WSS1_1DMaxGainPlusLast(this, a, dW, n, i)
            % Use last one as second index (algorithm can not happen to choose same indices for
            % two successive optimizations, as the analytical max is gained already.)
            in = this.WSS0_1D_GetMaxGainIndex(a, dW, n);
            if i ~= -1
                j = i;
                i = in;
            else
                [dWs, sortidx] = sort(abs(dW));
                pick = 1;
                while sortidx(pick) == in
                    pick = pick + 1;
                end
                i = in;
                j = sortidx(pick);
            end
        end
        
        function [i, j] = WSS2_TwoSetMaxGain(this, a, dW, n)
            len = round(n/2);
            set1 = 1:len;
            set2 = set1(end)+1:n;
            i = this.WSS0_1D_GetMaxGainIndex(a([set1 set1+n]), dW([set1 set1+n]), len);
            if i > len
                i = i+len-1;
            end
            j = this.WSS0_1D_GetMaxGainIndex(a([set2 set2+n]), dW([set2 set2+n]), numel(set2));
            if j > numel(set2)
                j = j + numel(set2)-1;
            end
            j=j+len;
        end
        
        function [i, j] = WSS4_TwoMaxGradients(this, a, dW)
            tmp = max(dW,-dW);
            [M, i] = max(tmp);
            [M, j] = max(M-tmp);
        end
        
        function [a, dW, T] = I0_ColdStart(this, fxi)
            T = 0;
            a = zeros(1,2*length(fxi)); %1..n = a^+, n+1..2n=a^-
            dW = [fxi -fxi] - this.Eps;
        end
        
        function [ai, dW, T] = I1_ColdStartKernelRule(this, fxi)
            n = length(fxi);
            %ai(1:n) = 1; ai(n+1:2*n) = 0;
            ai = max(0,min(this.C,[fxi.*(fxi>=0) -fxi.*(fxi<0)]));
            ai = ai * 10 / n;
            a = ai(1:n)-ai(n+1:end);
            hlp = -a*this.K + fxi;
            dW = [hlp -hlp] - this.Eps;
            %T = a*this.K*a' + this.Eps*sum(ai) - fxi*a';
            T = -sum(ai .* dW);
        end       
    end
    
    methods(Static)
        function res = test_ScalarEpsSVR_SMO(version)
            % Performs a test of this class
            
            if nargin == 0
                version = 1;
            end
            
            x = -5:.1:5;
            fx = sinc(x)+.2*x;
%             fx = sinc(x);
            
            svr = general.regression.ScalarEpsSVR_SMO;
            svr.Version = version;
            %svr.Eps = 0.073648;
            svr.Eps = .05;
            svr.Lambda = 1/20;%1/20; % i.e. C=10 as in ScalarEpsSVR
            
            %kernel = kernels.PolyKernel(7);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(.5);
            svr.K = kernel.evaluate(x,x);

            [ai, svidx] = svr.computeKernelCoefficients(fx);
            sv = x(:,svidx);
            svfun = @(x)ai'*(kernel.evaluate(x,sv)');
            
            fsvr = svfun(x);
            
            fdiff = abs(fsvr(svidx)-fx(svidx));
            errors = find(fdiff > 1.01*svr.Eps);
            res = isempty(errors);
            
            % Plot approximated function
            figure;
            plot(x,fx,'r',x,[fx-svr.Eps; fx+svr.Eps],'r--');
            hold on;
            plot(x,fsvr,'b',x,[fsvr-svr.Eps; fsvr+svr.Eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');
            
            if ~res
                plot(x(svidx(errors)),fx(svidx(errors)),'blackx','LineWidth',4);
            end
            
            tit = sprintf('#SV=%d, eps=%f',length(svidx),svr.Eps);
            title(tit);
            disp(tit);
            hold off;
        end
    end
end