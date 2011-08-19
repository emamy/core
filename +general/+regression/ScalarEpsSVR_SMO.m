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
        
        Lambda = 1;
        
        Vis = 1;
        
        Version = 2;
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
            C = 1/(2*this.Lambda);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            Sall = zeros(1,maxcnt);
            S = n*C;
            
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
                
                %% Get max gain
                % get clipped new \alpha^{+,-}
%                 ainew = max(0,min(C, dW + ai));%sgn.*
                
                % Check which alpha values might be changed (only the ones with their partner
                % variable equal to zero)
                ind = 1:n;
                ch = [ind(a(n+1:end) == 0) ind(a(1:n) == 0)+n];
                
                %% Find optima for changeable alphas
                r = max(-a(ch),min(C-a(ch),dW(ch)));%sgn(ch).* +this.Eps
                
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
                a(idx) = a(idx) + r;                               
                
                %% Get updated stopping criteria
                T = T + r*(r - 2*dW(idx) + sgn(idx)*fxi(idx1) - this.Eps);

                E = C * sum(max(0,min(2-this.Eps,dW)));
                
                %% Update dW - only on side with change
                Ki = this.K(idx1, :);
                if idx <= n
%                     dW(1:n) = dW(1:n) - r * Ki;
                    dW = dW - r * [Ki -Ki];
                else
%                     dW(n+1:end) = dW(n+1:end) - r * Ki;
                    dW = dW - r * [-Ki Ki];
                end
                
                %% Duality gap P-W (testing)
                dif = (a(1:n)-a(n+1:end)) * this.K;
                hlp = abs(fxi - dif) - this.Eps;
                hlp(hlp < 0) = 0;
                t_tmp = (dif-fxi)*(a(1:n)-a(n+1:end))' + this.Eps*sum(a);
                e_tmp = C * sum(hlp);
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
            this.Vis = 2;
            maxcnt = 1000;
            cnt = 1;
            
            n = length(fxi);
            C = 1/(2*this.Lambda);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            Sall = zeros(1,maxcnt);
            S = n*C;

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
                
                if cnt == 6
                    tool = 45;
                end
                
                %% Max gain computation - O(n^2) strategy
                [i, j] = this.WSS0_BestGainIndices(a, dW, C, sgn);
                %[i, j] = this.WSS1024_RandomAdmissibleIndices(a);
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
                r = (dW(i) - sgn(i) * kij * dW(j)) / (1-kij^2);
                s = (dW(j) - sgn(j) * kij * dW(i)) / (1-kij^2);
                
                % Handle constraint cases
                r2 = []; s2 = [];
                if r > C-a(i) || r < -a(i)
                    r = min(C-a(i), max(-a(i), r));
                    s2 = min(C-a(j), max(-a(j), dW(j) - r * kij));
                end
                if s > C-a(j) || s < -a(j)
                    s = min(C-a(j), max(-a(j), s));
                    r2 = min(C-a(i), max(-a(i), dW(i) - s * kij));
                end
                % Both updates are constrained, search for max gain
                if ~isempty(r2) && ~isempty(s2)
                    % Gain if r is constrained and s accordingly updated
                    gr = r*(dW(i)-.5*r) + s2*(dW(j)-.5*s2) - r*s2*kij;
                    % Gain if s is constrained and r accordingly updated
                    gs = r2*(dW(i)-.5*r2) + s*(dW(j)-.5*s) - r2*s*kij;
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
                
                if any(a < 0) || any(a > C)
                    warning('some:id','constraint violation. please check.');
                    keyboard;
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
                
                %% update stopping crits
                T = T + r*(r - 2*dW(i) - this.Eps + sgn(i)*fxi(i1)) ...
                      + s*(s - 2*dW(j) - this.Eps + sgn(j)*fxi(j1));

                E = C * sum(max(0,min(2-this.Eps,dW)));
                
                %% Duality gap P-W (testing)
                dif = (a(1:n)-a(n+1:end)) * this.K;
                hlp = abs(fxi - dif) - this.Eps;
                hlp(hlp < 0) = 0;
                t_tmp = (dif-fxi)*(a(1:n)-a(n+1:end))' + this.Eps*sum(a);
                e_tmp = C * sum(hlp);
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
        
        function [i, j] = WSS0_BestGainIndices(this, a, dW, C, sgn)
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
            ub = r > C - ai;
            lb = r < - ai;
                        
            % r is constrained only
            ronly = (ub | lb) & ~(ub' | lb');
            r(ronly) = min(C - ai(ronly),max(- ai(ronly), r(ronly)));
            s(ronly) = min(C-aj(ronly), max(-aj(ronly), wj(ronly) - r(ronly) .* Kch(ronly)));
            
            % s is constrained only
            sonly = (ub' | lb') & ~(ub | lb);
            s(sonly) = min(C - aj(sonly),max(- aj(sonly), s(sonly)));
            r(sonly) = min(C-ai(sonly), max(-ai(sonly), wi(sonly) - s(sonly) .* Kch(sonly)));
            
            % Case of both vars constrained
            % Update s2 and r2 for the cases of the other one is constrained
            both = (ub | lb) & (ub' | lb');
            r1 = zeros(size(r)); r2=r1; s1=r1; s2=r1;
            r1(both) = min(C - ai(both),max(- ai(both), r(both)));
            s1(both) = min(C-aj(both), max(-aj(both), wj(both) - r1(both) .* Kch(both)));
            s2(both) = min(C - aj(both),max(- aj(both), s(both)));
            r2(both) = min(C-ai(both), max(-ai(both), wi(both) - s2(both) .* Kch(both)));
            
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
        
        function [a, dW, T] = I0_ColdStart(this, fxi)
            T = 0;
            a = zeros(1,2*length(fxi)); %1..n = a^+, n+1..2n=a^-
            dW = [fxi -fxi] - this.Eps;
        end
        
        function [ai, dW, T] = I1_ColdStartKernelRule(this, fxi)
            n = length(fxi);
            C = 1/(2*this.Lambda);
            %ai(1:n) = 1; ai(n+1:2*n) = 0;
            ai = max(0,min(C,[fxi.*(fxi>=0) -fxi.*(fxi<0)]));
            ai = ai * 10 / n;
            a = ai(1:n)-ai(n+1:end);
            hlp = -a*this.K + fxi;
            dW = [hlp -hlp] - this.Eps;
            %T = a*this.K*a' + this.Eps*sum(ai) - fxi*a';
            T = -sum(ai .* dW);
        end
        
        function ai = regress1D_exp(this, fxi)
            % experimental version with direct negative update possibility
            
            %% Preps
            this.Vis = 0;
            this.Vis = 1;
            this.Vis = 2;
            maxcnt = 5000;
            cnt = 1;
            
            n = length(fxi);
            C = 1/(2*this.Lambda);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            S = zeros(1,maxcnt);
            S(cnt) = n*C;
            
            %% Init - cold start
            [a, dW, T] = I0_ColdStart(this, fxi);
            
            %% Init - kernel rule
            [a, dW, T] = I1_ColdStartKernelRule(this, fxi);
            
            if this.Vis > 0
                if this.Vis > 1
                    h = figure(1);
                end
                err = []; Tall = S; Tall(1) = T; Eall = S; Eall(1) = 0;
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
                r = max(-C-ai(ch)+2*this.Eps,min(C-ai(ch),dW(ch)));%sgn(ch).* +this.Eps
                
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
                
                if ai(idx) + r < 0
                   % switch sides case 
                   s = -(ai(idx) + r) - 2*this.Eps;%-2*this.Eps;
                   r = -ai(idx);
                   % switch from alpha^+ to alpha^-
                   if idx <= n
                       ai(idx+n) = ai(idx+n) + s;
                   % switch from alpha^- to alpha^+
                   else
                       ai(idx1) = ai(idx1) + s;
                   end
                else
                    s = 0;
                end
                
                %% update ai at max gain index
                ai(idx) = ai(idx) + r;
                
                %% Update dW - only on side with change
                Ki = this.K(idx1, :);
                if idx <= n
%                     dW(1:n) = dW(1:n) - r * Ki;
                    dW = dW - (r-s) * [Ki -Ki];
                else
%                     dW(n+1:end) = dW(n+1:end) - r * Ki;
                    dW = dW - (r-s) * [-Ki Ki];
                end
                
                %% update stopping crits
                yi = sgn(idx)*fxi(idx1);
                T = T + r*(r - 2*dW(idx) - this.Eps + yi);
                sidx = (idx <= n)*(idx+n) + (idx > n)*idx1;
                T = T + s*(s - 2*dW(sidx) - this.Eps + yi);
                
                dif = (ai(1:n)-ai(n+1:end)) * this.K;
                hlp = abs(fxi - min(1,max(-1,dif))) - this.Eps;
                hlp(hlp < 0) = 0;
                E = C * sum(hlp);

                % Unclipped E term
%                 hlp = abs(fxi - (ai(1:n)-ai(n+1:end)) * this.K) - this.Eps;
%                 hlp(hlp < 0) = 0;
%                 Eall = C * sum(hlp);

                cnt = cnt+1;
                S(cnt) = T + E;
                Tall(cnt) = T;
                Eall(cnt) = E;
            end
            
            fprintf('Finished after %d/%d iterations.\n',cnt,maxcnt);
            if this.Vis > 0
                figure;
                semilogy(1:cnt,S(1:cnt),'r',1:cnt,Tall(1:cnt),'b',1:cnt,Eall(1:cnt),'g',1:cnt,ones(1,cnt)*stop,'black');
                title('S, T, E and stopping values'); legend('S = T+E','T','E','stop cond');

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
% %                 Eall = C * sum(hlp)
% 
%                 cnt = cnt+1;
%                 S(cnt) = T + E;
%                 Tall(cnt) = T;
%                 Eall(cnt) = E;
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