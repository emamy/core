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
        
        StopEps = 1e-2;
        
        Lambda = 1;
        
        Vis = 1;
        
        Version = 1;
    end
    
    methods
        function ai = regress(this, fxi)
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
            maxcnt = 5000;
            cnt = 1;
            
            n = length(fxi);
            C = 1/(2*this.Lambda);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            S = zeros(1,maxcnt);
            S(cnt) = n*C;
            %% Init - cold start
            T = 0;
            ai = zeros(1,2*n); %1..n = a^+, n+1..2n=a^-
            dW = [fxi -fxi] - this.Eps;
            % End Init
            
            %% Init - kernel rule
            %1..n = a^+, n+1..2n=a^-
%             ai(1:n) = 1; ai(n+1:2*n) = 0;
%             %ai = max(0,min(C,[fxi.*(fxi>=0) -fxi.*(fxi<0)]))/n; 
%             a = ai(1:n)-ai(n+1:end);
%             T = a*this.K*a' + this.Eps*sum(ai) + fxi*a';
%             dW = [-ai(1:n)*this.K ai(n+1:end)*this.K] + [fxi -fxi] - this.Eps;
            % End Init
            
            if this.Vis > 0
                if this.Vis > 1
                    h = figure(1);
                end
                err = []; T2 = S; T2(1) = T; E2 = S; E2(1) = 0;
            end
            
            %rp = zeros(1:n);
            %rm = zeros(1:n);
            %dWp = dW(1:n);
            %dWm = dW(n+1:end);
            %ap = ai(1:n); am = ai(n+1:end);
            %g = zeros(1,1:2*n);
            
            stop = this.StopEps/(2*this.Lambda);
            while S(cnt) > stop && cnt < maxcnt
                
                if this.Vis > 0
                    afxi = (ai(1:n)-ai(n+1:end))*this.K;
                    if this.Vis > 1
                        figure(h);
                        subplot(1,2,1);
                        plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b');
                        title(sprintf('Current approximation, error: %e',err(end)));
                        legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                        axis tight;
                    end
                    etmp = abs(fxi-afxi);
                    err(end+1) = sum(etmp .* (etmp > this.Eps));%#ok
                    
                    %fprintf('%f ',sum(ai(1:n).*ai(n+1:end)));
                    if sum(ai(1:n).*ai(n+1:end)) ~= 0
                        warning('some:id','Invalid ai config. please check now!');
                        keyboard;
                    end
                end
                
                %% Get max gain
                % get clipped new \alpha^{+,-}
%                 ainew = max(0,min(C, dW + ai));%sgn.*
                
                % Check which alpha values might be changed (only the ones with their partner
                % variable equal to zero)
                ind = 1:n;
                ch = [ind(ai(n+1:end) == 0) ind(ai(1:n) == 0)+n];
                
                %% Find optima for changeable alphas
                r = max(-ai(ch),min(C-ai(ch),dW(ch)));%sgn(ch).* +this.Eps
                
                %% Working set selection
                % S1: Max W(\alpha) gain
                g = r .* (dW(ch) - .5*r);
                
                [M, idxch] = max(g);
                idx = ch(idxch);
                
                % S2: Max gradient value for nonbounded \alpha
%                 g = zeros(1,2*n);
%                 [M, idx] = max(max([dW .* (ai < C); -dW .* (ai > 0)],[],1));
                
                %% Misc vars
                % For index: Subtract n if index is for an \alpha^- (linear indexing)
                idx1 = idx - n*(idx > n);
                
                if this.Vis > 1
                    subplot(1,2,2);
                    r2 = zeros(1,2*n);
                    r2(ch) = r;
                    plot(r2,'b');
                    hold on;
                    g2 = zeros(1,2*n);
                    g2(ch) = g;
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
                
                %% select delta update for max gain index
                r = r(idxch);
                
                %% update ai at max gain index
                ai(idx) = ai(idx) + r;                               
                
                %% Update dW - only on side with change
                Ki = this.K(idx1, :);
                if idx <= n
%                     dW(1:n) = dW(1:n) - r * Ki;
                    dW = dW - r * [Ki -Ki];
                else
%                     dW(n+1:end) = dW(n+1:end) - r * Ki;
                    dW = dW - r * [-Ki Ki];
                end
                
                %% update stopping crits
                yi = sgn(idx)*fxi(idx1);
                T = T + r*(r - 2*dW(idx) - this.Eps + yi);
                
                dif = (ai(1:n)-ai(n+1:end)) * this.K;
                hlp = abs(fxi - min(1,max(-1,dif))) - this.Eps;
                hlp(hlp < 0) = 0;
                E = C * sum(hlp);

                % Unclipped E term
%                 hlp = abs(fxi - (ai(1:n)-ai(n+1:end)) * this.K) - this.Eps;
%                 hlp(hlp < 0) = 0;
%                 E2 = C * sum(hlp);

                cnt = cnt+1;
                S(cnt) = T + E;
                T2(cnt) = T;
                E2(cnt) = E;
            end
            
            fprintf('Finished after %d/%d iterations.\n',cnt,maxcnt);
            if this.Vis > 0
                figure;
                semilogy(1:cnt,S(1:cnt),'r',1:cnt,T2(1:cnt),'b',1:cnt,E2(1:cnt),'g',1:cnt,ones(1,cnt)*stop,'black');
                title('S, T, E and stopping values'); legend('S = T+E','T','E','stop cond');

                figure;
                semilogy(err);
                title('Approximation error (by eps-insensitive loss fcn)');
            end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (ai(1:n)-ai(n+1:end))';
        end
        
        function ai = regress2D(this, fxi)
            
            %% Preps
            this.Vis = 0;
            this.Vis = 1;
%             this.Vis = 2;
            maxcnt = 5000;
            cnt = 1;
            
            n = length(fxi);
            C = 1/(2*this.Lambda);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            S = zeros(1,maxcnt);
            S(cnt) = n*C;
            %% Init - cold start
            T = 0;
            ai = zeros(1,2*n); %1..n = a^+, n+1..2n=a^-
            dW = [fxi -fxi] - this.Eps;
            % End Init
            
            %% Init - kernel rule
            %1..n = a^+, n+1..2n=a^-
%             ai(1:n) = 1; ai(n+1:2*n) = 0;
%             %ai = max(0,min(C,[fxi.*(fxi>=0) -fxi.*(fxi<0)]))/n; 
%             a = ai(1:n)-ai(n+1:end);
%             T = a*this.K*a' + this.Eps*sum(ai) + fxi*a';
%             dW = [-ai(1:n)*this.K ai(n+1:end)*this.K] + [fxi -fxi] - this.Eps;
            % End Init
            
            if this.Vis > 0
                if this.Vis > 1
                    h = figure(1);
                end
                err = []; T2 = S; T2(1) = T; E2 = S; E2(1) = 0;
            end

            a = ai;

            stop = this.StopEps/(2*this.Lambda);
            while S(cnt) > stop && cnt < maxcnt
                
                if this.Vis > 0
                    afxi = (a(1:n)-a(n+1:end))*this.K;
                    etmp = abs(fxi-afxi);
                    err(end+1) = sum(etmp .* (etmp > this.Eps));%#ok
                    if this.Vis > 1
                        figure(h);
                        subplot(1,2,1);
                        plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b');
                        title(sprintf('Current approximation, error: %e',err(end)));
                        legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                        axis tight;
                    end
                    
                    %fprintf('%f ',sum(ai(1:n).*ai(n+1:end)));
                    if sum(a(1:n) .* a(n+1:end)) ~= 0
                        warning('some:id','Invalid ai config. please check now!');
                        keyboard;
                    end
                end
                
                %% Max gain computation - O(n^2) strategy
                % Checks which alpha values might be changed (only the ones with their partner
                % variable equal to zero)
                ind = 1:n;
                
                %% Find optima for changeable alpha^{+,-}
                apch = ind(a(n+1:end) == 0);                
                amch = ind(a(1:n) == 0);
                ch = [apch amch];
                Kch = this.K(ch,ch);
                ch(numel(apch)+1:end) = ch(numel(apch)+1:end)+n;
                [wj, wi] = meshgrid(dW(ch));
                ai = repmat(a(ch), length(ch),1);
                aj = repmat(a(ch)',1, length(ch));
                sig = sgn(ch)'*sgn(ch);
                div = 1./(1-Kch.^2);
                div(isinf(div)) = 0;
                r = (wi - sig .* Kch .* wj) .* div ;
                s = (wj - sig .* Kch .* wi) .* div;
%                 %% Kill NaN entries
%                 m = numel(ch); tru = true(m,1);
%                 setzero = spdiags([tru tru tru],[-numel(apch)+1, 0, numel(amch)],m,m);
%                 r(setzero) = 0;
%                 s(setzero) = 0;
                %% Handle constraints
                rup = r > C - ai;
                rlb = r < - ai;
                sup = s > C - aj;
                slb = s < - aj;
                % r is constrained
                r(rup) = C - ai(rup); r(rlb) = -ai(rlb);
                hlp = rup | rlb;
                s(hlp) = min(C-aj(hlp), max(-aj(hlp), wj(hlp) - r(hlp) .* Kch(hlp)));
                r2 = r; s2 = s;
                % s is constrained
                s(sup) = C - aj(sup); s(slb) = - aj(slb);
                hlp = sup | slb;
                r(hlp) = min(C-ai(hlp), max(-ai(hlp), wi(hlp) - s(hlp) .* Kch(hlp)));
                % r and s are both bounded
                % Determine update via get higher gain
                hlp = (rup | rlb) & (sup | rlb);
                g1 = r(hlp).*(wi(hlp)-.5*r(hlp)) + s(hlp).*(wj(hlp)-.5*s(hlp)) - r(hlp).*s(hlp).*Kch(hlp);
                g2 = r2(hlp).*(wi(hlp)-.5*r2(hlp)) + s2(hlp).*(wj(hlp)-.5*s2(hlp)) - r2(hlp).*s2(hlp).*Kch(hlp);
                cmp = g1 < g2;
                r(cmp) = r2(cmp); s(cmp) = s2(cmp);
                % Compute gain for r_i,s_j combinations
                g = r .* (wi - .5 * r) + s .* (wj - .5 * s) - r .* s .* Kch;
                % Extract i,j indices
                [hlp, jch] = max(g,[],1);
                [maxg, ich] = max(hlp);
                jch = jch(ich);
                i = ch(ich);
                j = ch(jch);
                r = r(ich, jch);
                s = s(ich, jch);

                a(i) = a(i) + r;
                a(j) = a(j) + s;
                
                if any(a) < 0 || any(a) > C
                    warning('some:id','constraint violation. please check.');
                    keyboard;
                end
             
                i1 = i - (i > n)*n;
                j1 = j - (j > n)*n;
                               
                if this.Vis > 1
                    [X,Y] = meshgrid(1:2*n,1:2*n);
                    subplot(1,2,2);
                    G = zeros(2*n,2*n);
                    G(ch,ch) = g;
                    surf(X,Y,G,'EdgeColor','none');
                    hold on;
                    plot3(X(ch(ich),ch(jch)),Y(ch(ich),ch(jch)),g(ich,jch),'r.','MarkerSize',6);
                    hold off;
                    title(sprintf('Best gain index: (%d,%d), gain: %e, \\alpha_{%d} change: %e, \\alpha_{%d} change: %e',i,j, maxg, i, r, j,s));
                    axis tight;
                    subplot(1,2,1);
                    hold on;
                    plot([i1 j1],afxi([i1 j1]),'.','MarkerSize',5);
                    hold off;
                end
                
                %% Update dW - only on side with change
                Ki = this.K(i1, :);
                Kj = this.K(j1, :);
                
                % Update dW
                if i <= n
                    dW = dW - r * [Ki -Ki];
                else
                    dW = dW - r * [-Ki Ki];
                end
                if j <= n
                    dW = dW - s * [Kj -Kj];
                else
                    dW = dW - s * [-Kj Kj];
                end

%                 if vp >= vm
%                     dWp = dWp - r * Ki - s * Kj;
%                     dWm = dWm + r * Ki + s * Kj;
%                     dW = dWp;
%                 else
%                     dWm = dWm - r * Ki - s * Kj;
%                     dWp = dWp + r * Ki + s * Kj;
%                     dW = dWm;
%                 end
                
                %% update stopping crits
                T = T + r*(r - 2*dW(i) - this.Eps + sgn(i)*fxi(i1)) ...
                      + s*(s - 2*dW(j) - this.Eps + sgn(j)*fxi(j1));
                
                dif = (a(1:n)-a(n+1:end)) * this.K;
                hlp = abs(fxi - min(1,max(-1,dif))) - this.Eps;
                hlp(hlp < 0) = 0;
                E = C * sum(hlp);

                % Unclipped E term
%                 hlp = abs(fxi - (ai(1:n)-ai(n+1:end)) * this.K) - this.Eps;
%                 hlp(hlp < 0) = 0;
%                 E2 = C * sum(hlp);

                cnt = cnt+1;
                S(cnt) = T + E;
                T2(cnt) = T;
                E2(cnt) = E;
%                 pause;
            end
            
            fprintf('Finished after %d/%d iterations.\n',cnt,maxcnt);
            if this.Vis > 0
                figure;
                semilogy(1:cnt,S(1:cnt),'r',1:cnt,T2(1:cnt),'b',1:cnt,E2(1:cnt),'g',1:cnt,ones(1,cnt)*stop,'black');
                title('S, T, E and stopping values'); legend('S = T+E','T','E','stop cond');

                figure;
                semilogy(err);
                title('Approximation error (by eps-insensitive loss fcn)');
            end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (a(1:n)-a(n+1:end))';
        end
        
        function ai = regress1D_testing(this, fxi)
            
            %% Preps
%             showVis = false;
            showVis = true;
            maxcnt = 5000;
            cnt = 1;
            
            n = length(fxi);
            C = 1/(2*this.Lambda);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            S = zeros(1,maxcnt);
            S(cnt) = n*C;
            %% Init - cold start
            T = 0;
            ai = zeros(1,n); %1..n = a^+, n+1..2n=a^-
            dW = [fxi -fxi] - this.Eps;
            % End Init
            
            %% Init - kernel rule
            %1..n = a^+, n+1..2n=a^-
%             ai(1:n) = 1; ai(n+1:2*n) = 0;
%             %ai = max(0,min(C,[fxi.*(fxi>=0) -fxi.*(fxi<0)]))/n; 
%             a = ai(1:n)-ai(n+1:end);
%             T = a*this.K*a' + this.Eps*sum(ai) + fxi*a';
%             dW = [-ai(1:n)*this.K ai(n+1:end)*this.K] + [fxi -fxi] - this.Eps;
            % End Init
            
            if showVis
                h = figure(1);
                err = []; T2 = S; T2(1) = T; E2 = S; E2(1) = 0;
            end
            
            %rp = zeros(1:n);
            %rm = zeros(1:n);
            %dWp = dW(1:n);
            %dWm = dW(n+1:end);
            %ap = ai(1:n); am = ai(n+1:end);
            %g = zeros(1,1:2*n);
            
            stop = this.StopEps/(2*this.Lambda);
            while S(cnt) > stop && cnt < maxcnt
                
                if showVis
                    figure(h);
                    subplot(1,2,1);
                    afxi = (ai(1:n)-ai(n+1:end))*this.K;
                    plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b');
                    etmp = abs(fxi-afxi);
                    err(end+1) = sum(etmp .* (etmp > this.Eps));
                    title(sprintf('Current approximation, error: %e',err(end)));
                    legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                    axis tight;
                    
                    %fprintf('%f ',sum(ai(1:n).*ai(n+1:end)));
                    if sum(ai(1:n).*ai(n+1:end)) ~= 0
                        warning('some:id','Invalid ai config. please check now!');
                        keyboard;
                    end
                end
                
                %% Get max gain
                % get clipped new \alpha^{+,-}
%                 ainew = max(0,min(C, dW + ai));%sgn.*
                
                % Check which alpha values might be changed (only the ones with their partner
                % variable equal to zero)
%                 ind = 1:n;
%                 ch = [ind(ai(n+1:end) == 0) ind(ai(1:n) == 0)+n];
                
                %% Find optima for changeable alphas
                r = max(-ai(ch),min(C-ai(ch),dW(ch)+this.Eps));%sgn(ch).* +this.Eps
                
                %% Working set selection
                % S1: Max W(\alpha) gain
                g = r .* (dW(ch) - .5*r);
                
                [M, idxch] = max(g);
                idx = ch(idxch);
                % S2: Max gradient value for nonbounded \alpha
%                 g = zeros(1,2*n);
%                 [M, idx] = max(max([dW .* (ai < C); -dW .* (ai > 0)],[],1));
                
                %% Misc vars
                % For index: Subtract n if index is for an \alpha^- (linear indexing)
                idx1 = idx - n*(idx > n);
                
                if showVis
                    subplot(1,2,2);
                    r2 = zeros(1,2*n);
                    r2(ch) = r;
                    plot(r2,'b');
                    hold on;
                    g2 = zeros(1,2*n);
                    g2(ch) = g;
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
                
                %% select delta update for max gain index
                r = r(idxch);
                
                %% update ai at max gain index
                ai(idx) = ai(idx) + r;                               
                
                %% Update dW - only on side with change
                Ki = this.K(idx1, :);
                if idx <= n
                    dW(1:n) = dW(1:n) - r * Ki;
                else
                    dW(n+1:end) = dW(n+1:end) - r * Ki;
                end
                
                %% update stopping crits
                yi = sgn(idx)*fxi(idx1);
                T = T + r*(r - 2*dW(idx) - this.Eps + yi);
                
                hlp = abs(fxi - min(1,max(-1,(ai(1:n)-ai(n+1:end)) * this.K))) - this.Eps;
                hlp(hlp < 0) = 0;
                E = C * sum(hlp);

                % Unclipped E term
%                 hlp = abs(fxi - (ai(1:n)-ai(n+1:end)) * this.K) - this.Eps;
%                 hlp(hlp < 0) = 0;
%                 E2 = C * sum(hlp)

                cnt = cnt+1;
                S(cnt) = T + E;
                T2(cnt) = T;
                E2(cnt) = E;
            end
            
            fprintf('Finished after %d/%d iterations.\n',cnt,maxcnt);
            if showVis
                figure;
                semilogy(1:cnt,S(1:cnt),'r',1:cnt,T2(1:cnt),'b',1:cnt,E2(1:cnt),'g');
                title('S, T, E values'); legend('S = T+E','T','E');

                figure;
                semilogy(err);
                title('Approximation error (by eps-insensitive loss fcn)');
            end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (ai(1:n)-ai(n+1:end))';
        end
    end
    
    methods(Static)
        function res = test_ScalarEpsSVR_SMO
            % Performs a test of this class
            
            x = -5:.1:5;
            fx = sinc(x)+.2*x;
%             fx = sinc(x);
            
            svr = general.regression.ScalarEpsSVR_SMO;
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

 %                 % tmp test
    %                 % DIRECT ainew computation NOT using gradients (fits, put inside while loop)
    %                 ainew2 = zeros(size(ainew));
    %                 for p = 1:n
    %                     ap = ai(1:n); ap(p) = 0;
    %                     ainew2(p) = fxi(p) - this.Eps - (ap-ai(n+1:end)) * this.K(:,p);
    %                     am = ai(n+1:end); am(p) = 0;
    %                     ainew2(n+p) =  - fxi(p) - this.Eps + (ai(1:n)-am) * this.K(:,p);
    %                 end
    %                 ainew2 = max(0,min(C, ainew2));
    %                 d2 = ainew2 - ai;
    %                 g2 = d2 .* (dW - .5*d2);
    %                 % tmp test end

%         while S(cnt) > stop && cnt < maxcnt
%                 %% Get max gain
%                 % get clipped new \alpha^{+,-}
% %                 ainew = max(0,min(C, dW + ai));%sgn.*
% 
%                 % Test ainew comp for ai^+
%                 ap = ai(1:n); am = ai(n+1:end);
%                 delta = fxi - this.Eps - (ap-am)*this.K;
%                 r1 = -.5*(ap+am-delta) + .5 * sign(delta) .* sqrt(am.^2 + ap.^2 + delta.^2 + 2*delta .* (ap-am));
%                 %r1b = -.5*(ap+am-delta) - .5*sqrt(am.^2 + ap.^2 + delta.^2 + 2*delta .* (ap-am));
%                 ainew(1:n) = ap + r1;
%                 %ainew2(1:n) = ap + r2;
%                 dp1 = (r1.*am)./(ap + r1); 
%                 dp1(isnan(dp1)) = 0;
%                 %dp2 = r1b + (r1b.*am)/(ap + r1b); 
%                 
%                 delta2 = (ap-am)*this.K - this.Eps - fxi;
%                 r2 = -.5*(ap+am-delta2) + .5 * sign(delta2) .* sqrt(am.^2 + ap.^2 + delta2.^2 - 2*delta2 .* (ap-am));
%                 %r2b = -.5*(ap+am-delta2) - .5*sqrt(am.^2 + ap.^2 + delta2.^2 - 2*delta2 .* (ap-am));
%                 
%                 ainew(n+1:2*n) = am + r2;
%                 %ainew2(n+1:2*n) = am + r22;
%                 dm1 = (r2.*ap)./(am + r2);
%                 dm1(isnan(dm1)) = 0;
%                 %dm2 = r2b + (r2b.*ap)/(am + r2b); 
%                 
%                 % get delta \alpha^{new} - \alpha
%                 %d = ainew - ai;
%                 d = [dp1+r1 dm1+r2];
%                 %d2 = [dp2 dm2];
%                 
%                 %% Working set selection
%                 % S1: Max W(\alpha) gain
%                 g = d .* (dW - .5*d);
%                 %g2 = d2 .* (dW - .5*d2);
%                 
%                 [M, idx] = max(g);
%                 %[M2, idx2] = max(g2);
%                 % S2: Max gradient value for nonbounded \alpha
% %                 g = zeros(1,2*n);
% %                 [M, idx] = max(max([dW .* (ai < C); -dW .* (ai > 0)],[],1));
%                 
%                 %% Misc vars
%                 % For index: Subtract n if index is for an \alpha^- (linear indexing)
%                 idx1 = idx - n*(idx > n);
%                 
%                 if showVis
%                     figure(h);
%                     subplot(1,2,1);
%                     plot(d(1:n)-d(n+1:end),'b');
%                     hold on;
%                     plot(sqrt(g(1:n))-sqrt(g(n+1:end)),'g');
%                     plot(idx1,sgn(idx)*sqrt(g(idx)),'r.','MarkerSize',6);
%                     hold off;
%                     title(sprintf('Best gain index: %d, gain: %e, \\alpha_{%d} change: %e',idx1, M, idx1, d(idx)));
%                     legend('\alpha difference','sqrt(gain)');
%                     axis tight;
%                     
%                     subplot(1,2,2);
%                     afxi = (ai(1:n)-ai(n+1:end))*this.K;
%                     plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b',idx1,afxi(idx1),'.','MarkerSize',5);
%                     etmp = abs(fxi-afxi);
%                     err(end+1) = sum(etmp .* (etmp > this.Eps));
%                     title(sprintf('Current approximation, error: %e',err(end)));
%                     legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
%                     axis tight;
% %                     pause;
%                 end
%                 
%                 %% select delta update for max gain index
%                 d = d(idx);
%                 
%                 %% update ai at max gain index
%                 ai(idx) = ainew(idx);
%                 if idx <= n
%                     ai(idx+n) = ai(idx+n) - dp1(idx);
%                 else
%                     ai(idx1) = ai(idx1) - dm1(idx1);
%                 end
%                 ai = min(C,max(0,ai));
%                 
%                 chksum = sum(ai(1:n) .* ai(n+1:end))
%                 
%                 %% Update overall dW
%                 Ki = this.K(idx1, :);
%                 dW = dW - d * [Ki Ki];
%                 
%                 %% update stopping crits
%                 yi = sgn(idx)*fxi(idx1);
%                 T = T + d*(d - 2*dW(idx) - this.Eps + yi);
%                 
%                 hlp = abs(fxi - min(1,max(-1,(ai(1:n)-ai(n+1:end)) * this.K))) - this.Eps;
%                 hlp(hlp < 0) = 0;
%                 E = C * sum(hlp);
% 
%                 % Unclipped E term
% %                 hlp = abs(fxi - (ai(1:n)-ai(n+1:end)) * this.K) - this.Eps;
% %                 hlp(hlp < 0) = 0;
% %                 E2 = C * sum(hlp)
% 
%                 cnt = cnt+1;
%                 S(cnt) = T + E;
%                 T2(cnt) = T;
%                 E2(cnt) = E;
%             end

% %% Handle upper constraints
%                 % r_i > C-a_i
%                 rup = r > C - aip;
%                 % s_j > C-a_j
%                 sup = s > C - ajp;
%                 % r is constrained above
%                 r_rup = zeros(size(r)); s_rup = r_rup;
%                 r_rup(rup) = C - aip(rup);
%                 s_rup(rup) = wj(rup) - r_rup(rup) .* Kch(rup);
%                 % s is constrained above
%                 r_sup = zeros(size(r)); s_sup = r_sup;
%                 s_sup(sup) = C-ajp(sup);
%                 r_sup(sup) = wi(sup) - s_sup(sup) .* Kch(sup);
%                 % Transfer updates
%                 % - r upper bounded (all, incl. s upper bounded)
%                 r(rup) = r_rup(rup); s(rup) = s_rup(rup);
%                 % - s upper bounded (all, incl. r upper bounded)
%                 r(sup) = r_sup(sup); s(sup) = s_sup(sup);
%                 % - r,s upper bounded (overwrites possibly wrong values from above cases)
%                 % Determine update via get higher gain
%                 hlp = rup & sup;
%                 g1 = r_rup(hlp).*(wi(hlp)-.5*r_rup(hlp)) + s_rup(hlp).*(wj(hlp)-.5*s_rup(hlp)) - r_rup(hlp).*s_rup(hlp).*Kch(hlp);
%                 g2 = r_sup(hlp).*(wi(hlp)-.5*r_sup(hlp)) + s_sup(hlp).*(wj(hlp)-.5*s_sup(hlp)) - r_sup(hlp).*s_sup(hlp).*Kch(hlp);
%                 cmp = g2 <= g1;
%                 r(cmp) = r_rup(cmp); s(cmp) = s_rup(cmp);
%                 r(~cmp) = r_sup(cmp); s(~cmp) = s_sup(cmp);