classdef ScalarEpsSVR_SMO < general.regression.BaseScalarSVR
% ScalarEpsSVR_SMO: 
%
% Implements the 1D and 2D SMO algorithm introduced in SHS11: 
% 'I. Steinwart, D. Hush, and C. Scovel. Training svms without o?set. J. Mach. Learn. Res.,
% 12:141?202, February 2011'
%
% @author Daniel Wirtz @date 2011-08-09
%
% @change{0,5,dw,2011-09-09} Finished initial work on this class, fixed many bugs and implemented
% all strategies required for the WSS7 described in SHS11. Added initial
% value (warm start) support.
%
% @new{0,5,dw,2011-08-09} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The epsilon value to use for the loss function.
        %
        % @propclass{critical} Influences the precision of the result.
        %
        % @default .1;
        Eps = .1;
        
        % Stopping criterion.
        %
        % @propclass{critical} Determines the precision of the solution. This eps is used within the
        % stopping criterion introduced in SHS11.
        StopEps = 1e-4;
    
        % Visual output.
        %
        % @propclass{optional}
        %
        % @default 0;
        Vis = 0;
        
        % Which Version to use.
        %
        % Set to 1 for 1D SMO and 2 for 2D SMO.
        % 2D is faster, uses the WSS7 as default.
        %
        % @propclass{important} Influences the performance and time needed for solving.
        %
        % @default 2
        Version = 2;
        
        % Nearest neighbors to search for WSS4 strategy.
        %
        % @propclass{optional} A higher number of neighbors might 
        %
        % @default 10
        NNk = 10;
        
        % The maximum number of iterations.
        %
        % @propclass{alglimit} Max iteration count.
        %
        % @default 60000
        MaxCount = 60000;
    end
    
    methods
        function ai = regress(this, fxi, initialai)
            % Make sure it's a row vector.
            fxi = reshape(fxi,1,[]);
%             initialai = reshape(initialai,1,[]);
            if this.Version == 1
                % Call 1D-Optimizer
                ai = this.regress1D(fxi, initialai);
            else
                % Call 2D-Optimizer
                ai = this.regress2D(fxi, initialai);
            end
        end
    end
    
    methods(Access=private)
        
        function ai = regress1D(this, fxi, initialai)
            
            if ~isempty(initialai)
                % Init - warm start
                [a, dW, T] = I2_WarmStart(this, fxi, initialai);
            else
                % Init - cold start
                [a, dW, T] = I0_ColdStart(this, fxi);
                
                % Init - kernel rule
%                 [a, dW, T] = I1_ColdStartKernelRule(this, fxi);
            end

            %% Preps
            cnt = 1;
            
            n = length(fxi);
            
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            S = n*this.C;
            
            if this.Vis > 1
                h = figure(1);
            end
            
            stop = this.StopEps/(2*this.Lambda);
            while S > stop && cnt < this.MaxCount
                
                if this.Vis > 1
                    afxi = (a(1:n)-a(n+1:end))*this.K;
                    figure(h);
                    subplot(1,2,1);
                    plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps],'r--',1:n,afxi,'b');
                    %err = sum((fxi-afxi).^2);
                    %title(sprintf('Current approximation, error: %e',err(cnt)));
                    legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                    axis tight;
                end
                
                %% Select working set
                [idx, M] = this.WSS0_1D_GetMaxGainIndex(a, dW);
                
                %[idx, M] = this.WSS1_1D_GetMaxGradientIndex(a, dW);
                
                r = max(-a(idx),min(this.C-a(idx),dW(idx)));
                
                %% Misc vars
                % For index: Subtract n if index is for an \alpha^- (linear indexing)
                idx1 = idx - n*(idx > n);
                
                if this.Vis > 1
                    subplot(1,2,1);
                    hold on;
                    plot(idx1,afxi(idx1),'.','MarkerSize',5);
                    hold off;
                end
                
                %% update ai at max gain index
                a(idx) = a(idx) + r;                  
                
                %% Update T with old dW
                T = T + r*(r - 2*dW(idx) + sgn(idx)*fxi(idx1) - this.Eps);

                %% Update dW - only on side with change
                Ki = this.K(idx1, :);
                dW = dW + sign(idx-n-.5) * r * [Ki -Ki];
                
                %% Update E
                E = this.C * sum(max(0,min(2-this.Eps,dW)));
                                
                S = T + E;
                cnt = cnt+1;
            end
            if this.Vis > 0
                fprintf('Finished after %d/%d iterations.\n',cnt,this.MaxCount);
            end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (a(1:n)-a(n+1:end))';
        end
        
        function ai = regress2D(this, fxi, initialai)
            
            if ~isempty(initialai)
                % Init - warm start
                [a, dW, T] = I2_WarmStart(this, fxi, initialai);
            else
                % Init - cold start
                [a, dW, T] = I0_ColdStart(this, fxi);

                % Init - kernel rule
%                 [a, dW, T] = I1_ColdStartKernelRule(this, fxi);
            end
            
            if this.Vis > 1
                h = figure(1);
            end
            
            %% Preps
            n = length(fxi);
            sgn = ones(1,2*n);
            sgn(n+1:end) = -1;
            
            S = n*this.C;

            % Initialize i to -1 if MaxGainPlusLast is used.
            i = -1;
            
            cnt = 1;
            stop = this.StopEps/(2*this.Lambda);
            while abs(S) > stop && cnt < this.MaxCount+1
                
                %% Max gain computation - O(n^2) strategy
                ij = [];%#ok<*AGROW>
                
%                 ij = [ij; this.WSS0_BestGainIndices(a, dW, sgn);
                
%                 ij = [ij; this.WSS1024_RandomAdmissibleIndices(a);
                
%                 ij = [ij; this.WSS1_1DMaxGainPlusLast(a, dW, i)];
                
%                 ij = [ij; this.WSS2_TwoSetMax1DGain(a, dW)];
                
%                 ij = [ij; this.WSS4_1DMaxGainPluskNN(a, dW)];
                
                ij = [ij; this.WSS7_Combi(a, dW, i)];
                
                [r, s, idx] = this.getMax2DGainUpdatesForIndices(a, dW, ij);
                
                % Not usable with current scheme
                %[i, j] = WSS128_ApproxBestGainIndices(this, a, dW, sgn);
                
                i = ij(idx,1);
                j = ij(idx,2);
                i1 = i - (i > n)*n;
                j1 = j - (j > n)*n;
                
                if this.Vis > 1
                    afxi = (a(1:n)-a(n+1:end))*this.K;
                    figure(h);
                    subplot(1,2,1);
                    plot(1:n,fxi,'r',1:n,[fxi-this.Eps; fxi+this.Eps]','r--',1:n,afxi,'b');
                    hold on;
                    if (i1 == i)
                        plot(i,afxi(i),'o','MarkerEdgeColor','k',...
                            'MarkerFaceColor',[.1 .8 .1],...
                            'MarkerSize',6);
                    else
                        plot(i1,afxi(i1),'o','MarkerEdgeColor','k',...
                            'MarkerFaceColor',[.8 .1 .1],...
                            'MarkerSize',6);
                    end
                    if (j1 == j)
                        plot(j,afxi(j),'o','MarkerEdgeColor','k',...
                            'MarkerFaceColor',[.1 .8 .1],...
                            'MarkerSize',6);
                    else
                        plot(j1,afxi(j1),'o','MarkerEdgeColor','k',...
                            'MarkerFaceColor',[.8 .1 .1],...
                            'MarkerSize',6);
                    end
                    
                    hold off;
%                     title(sprintf('Current approximation, error: %e',err(cnt)));
                    %legend('f(x_i) training values','+\epsilon','-\epsilon','approximation values');
                    axis tight;
                    if this.Vis > 2
                        fprintf('alpha_{%d} change: %e, alpha_{%d} change: %e\n',i,r,j,s);
                    end
                    pause;
                end
                
                a(i) = a(i) + r;
                a(j) = a(j) + s;
                
                % If difference is smaller than n, they are ++ or -- updates, meaning addition.
                % otherwise, is |i-j| \geq n, we have +- or -+ cases, meaning subtraction.
                % subtraction of .5 avoids having sign(0) = 0.
                pm = sign(i - n - .5)*sign(j - n - .5);
                %% Update T using old dW
                T = T + r*(r - 2*dW(i) - this.Eps + sgn(i)*fxi(i1)) ...
                      + s*(s - 2*dW(j) - this.Eps + sgn(j)*fxi(j1)) ...
                      + 2*pm*r*s*this.K(i1,j1);
                
                %% Update dW - only on side with change
                Ki = this.K(i1, :);
                Kj = this.K(j1, :);
                dW = dW + sign(i-n-.5) * r * [Ki -Ki] ...
                        + sign(j-n-.5) * s * [Kj -Kj];
                
                %% Get new E term
                E = this.C * sum(max(0,min(2-this.Eps,dW)));
                
                S = T + E;
%                 Sall(cnt) = S;
%                 Tall(cnt) = T;
%                 Eall(cnt) = E;
                cnt = cnt+1;
            end
            cnt = cnt-1;
            if this.Vis > 0 || KerMor.App.Verbose > 0
                fprintf('ScalarEpsSVR_SMO: Finished after %d/%d iterations.\n',cnt,this.MaxCount);
            end
%             if this.Vis > 0
%                 figure;
%                 semilogy(1:cnt,abs(Sall(1:cnt)),'r',1:cnt,abs(Tall(1:cnt)),'b',1:cnt,Eall(1:cnt),'g',1:cnt,ones(1,cnt)*stop,'black');
%                 axis tight;
%                 title('S, T, E and stopping values'); legend('S = T+E','T','E','stop cond');
% 
%                 figure;
%                 semilogy(err(1:cnt));
%                 axis tight;
%                 title('Approximation error (by eps-insensitive loss fcn)');
%             end
            
            % \alpha = \alpha^+ - \alpha^-
            ai = (a(1:n)-a(n+1:end))';            
        end
        
        function [i, M] = WSS0_1D_GetMaxGainIndex(this, a, dW)
            n = numel(a)/2;
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
            
            if this.Vis > 1
                subplot(1,2,2);
                r2 = zeros(size(a));
                r2(ch) = r;
                plot(1:2*n,r2,'b');
                hold on;
                g2 = zeros(size(a));
                g2(ch) = g;
                plot(1:2*n,g2,'g');
                plot(i,g(idxch),'r.','MarkerSize',6);
                plot([n n+eps],[0, max([r g])],'black');
                hold off;
                title(sprintf('Best gain index: %d, gain: %e, \\alpha_{%d} change: %e',i, M, i, r(idxch)));
                legend('\alpha difference','gain');
                axis tight;
            end
        end
        
        function [i, M] = WSS1_1D_GetMaxGradientIndex(this, a, dW)
             % S2: Max gradient value for nonbounded \alpha
             [M, i] = max(max([dW .* (a < this.C); -dW .* (a > 0)],[],1));
        end
        
        function ij = WSS0_AllIndices(this, a)%#ok
            % Implements the O(n^2)-WSS0 working set selection strategy.
            
            n = numel(a)/2;
            ind = 1:n;
            apch = ind(a(n+1:end) == 0);
            amch = ind(a(1:n) == 0);
            ch = [apch amch+n];
            % Create unique combinations of indices
            A = repmat(ch,length(ch),1);
            I = triu(A,1);
            J = triu(A',1);
            ij = [I(I~=0) J(J~=0)];
            % Remove indices that would look at moving the same alpha up and down at the same time
            ij(ij(:,1)-ij(:,2)-n == 0,:) = [];
        end
        
        function ij = WSS1_1DMaxGainPlusLast(this, a, dW, i)
            % Use last one as second index (algorithm can not happen to choose same indices for
            % two successive optimizations, as the analytical max is gained already.)
            in = this.WSS0_1D_GetMaxGainIndex(a, dW);
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
            ij = [i j];
        end
        
        function ij = WSS2_TwoSetMax1DGain(this, a, dW)
            %
            % Splitting:
            %         \alpha^+                \alpha^-
            %|-----------------------|-------------------------|
            %1                       n                         2n
            %    set1^+                   set1^-
            %|-----------|            |-----------|         
            %1           len          n+1         n+len
            %                set2^+                    set2^-
            %             |----------|             |-----------|         
            %1            len+1      n             n+len+1     2n
            n = numel(a)/2;
            len = round(n/2);
            set1 = 1:len;
            set2 = set1(end)+1:n;
            i = this.WSS0_1D_GetMaxGainIndex(a([set1 set1+n]), dW([set1 set1+n]));
            if i > len
                i = i+(n-len); % Move to \am index (add number of set2 elements), i.e. i = i+numel(set2)
            end
            j = this.WSS0_1D_GetMaxGainIndex(a([set2 set2+n]), dW([set2 set2+n]));
            if j > numel(set2)
                j = j+len; % If larger than set2 size add number of set1 elements
            end
            % Add +len as this is the offset for the second set in each case
            ij = [i j+len];
        end
        
        function ij = WSS4_1DMaxGainPluskNN(this, a, dW)
            i = this.WSS0_1D_GetMaxGainIndex(a, dW);
            
            n = numel(a)/2;
            ki = i;
            if ki > n
                ki = ki - n;
            end
            [dist, idx] = sort(sqrt(1-this.K(ki,:)));
            % Extract indices of alpha^+ that are "close" wrt to the kernel metric \sqrt{1-\Phi(x_i,x_j)}
            % Start at 2 as i itself is not an option for second index
            nidx = idx(2: min(this.NNk+1,numel(idx)));
            % Of course also consider the \alpha^-
            ij = [i*ones(size(nidx,2)*2,1) [nidx nidx+n]'];
        end
        
        function ij = WSS7_Combi(this, a, dW, i)
            ij = this.WSS2_TwoSetMax1DGain(a, dW);
                try
            ij = [ij; this.WSS4_1DMaxGainPluskNN(a, dW)];
                catch ME
                    keyboard;
                end
            
            % Implement Max1DGainPlusLast directly to save a call to the 1D method.
            in = ij(end,1);
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
            ij = [ij; i j];
        end
        
        function ij = WSS128_ApproxBestGainIndices(this, a, dW)
            % Implements the O(n^2)-WSS0 working set selection strategy.
            
            % Only consider changeable \alpha^{+,-}, i.e. the ones whos partner \alpha equals zero.
            n = numel(a)/2;
            ind = 1:n;
            apch = ind(a(n+1:end) == 0);
            amch = ind(a(1:n) == 0);
            ch = [apch amch+n];
            % Create unique combinations of indices
            A = repmat(ch,length(ch),1);
            I = triu(A,1);
            J = triu(A',1);
            ij = [I(I~=0) J(J~=0)];
            % Remove indices that would look at moving the same alpha up and down at the same time
            ij(ij(:,1)-ij(:,2)-n == 0,:) = [];
            
            kijidx = ij;
            kijidxup = kijidx > n;
            kijidx(kijidxup) = kijidx(kijidxup)-n;
            kijidx = sub2ind(size(this.K),kijidx(:,1),kijidx(:,2));
            Kij = this.K(kijidx);
            
            sgn = 1-mod(sum(kijidxup,2),2)*2;
            dwi = dW(ij);
            aij = a(ij);
            
            div = 1 ./ (1-Kij.^2);
            
            % Compute update r
            rs = [(dwi(:,1) - sgn .* Kij .* dwi(:,2)) .* div ...
                  (dwi(:,2) - sgn .* Kij .* dwi(:,1)) .* div];
            ub = rs > this.C - aij;
            lb = rs < - aij;
            rs(ub) = this.C-aij(ub);
            rs(lb) = -aij(lb);            
            
            %% Compute gain for r_i,s_j combinations
            g = sum(rs .* (dwi - .5*rs),2) - sgn.*rs(:,1).*rs(:,2).*Kij;
            
            % Extract i,j indices
            [maxg,idx] = max(g);
            
            if this.Vis > 2
                fprintf('Updates found inside WSS128: r=%f, s=%f\n',rs(idx,1),rs(idx,2));
            end
            ij = ij(idx,:); 
            if this.Vis > 1    
                subplot(1,2,2);
                plot(1:length(g),max(0,g));
                hold on;
                plot(idx,g(idx),'black.','MarkerSize',8);
                hold off;
                title(sprintf('Best gain index: (%d,%d), gain: %e',ij(1),ij(2), maxg));
                axis tight;
            end
        end
        
        function ij = WSS1024_RandomAdmissibleIndices(this, a)%#ok
            n = numel(a)/2;
            ind = 1:n;
            ch = [ind(a(n+1:end) == 0) ind(a(1:n) == 0)+n];
            i = ch(randi(length(ch)));
            j = ch(randi(length(ch)));
            while j==i || abs(j-i) == n
                j = ch(randi(length(ch)));
            end
            ij = [i j];
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
        
        function [ai, dW, T] = I2_WarmStart(this, fxi, initialai)
            n = length(fxi);
            ap = zeros(1,n);
            am = ap;
            ap(initialai > 0) = initialai(initialai>0);
            am(initialai < 0) = initialai(initialai<0);
            ai = [ap am];
            hlp = -(ap-am)*this.K + fxi;
            dW = [hlp -hlp] - this.Eps;
            T = -sum(ai .* dW);
        end
        
        function [r, s, idx] = getMax2DGainUpdatesForIndices(this, a, dW, ij)
            n = numel(a)/2;
            kijidx = ij;
            kijidxup = kijidx > n;
            kijidx(kijidxup) = kijidx(kijidxup)-n;
            kijidx = sub2ind([n n],kijidx(:,1),kijidx(:,2));
            Kij = this.K(kijidx);
            
            sgn = 1-mod(sum(kijidxup,2),2)*2;
            dwi = dW(ij);
            aij = a(ij);
            
            div = 1 ./ (1-Kij.^2);
            
            % Compute update r
            rs = [(dwi(:,1) - sgn .* Kij .* dwi(:,2)) .* div ...
                  (dwi(:,2) - sgn .* Kij .* dwi(:,1)) .* div];
            
            % Handle constraints
            b = rs > this.C - aij | rs < - aij;
            
            ronly = find(b(:,1) & ~b(:,2));
            sonly = find(~b(:,1) & b(:,2));
            both = find(b(:,1) & b(:,2));
            
            % r is constrained only
            if ~isempty(ronly)
                rs(ronly,1) = min(this.C - aij(ronly,1),max(- aij(ronly,1), rs(ronly,1)));
                rs(ronly,2) = min(this.C-aij(ronly,2), max(-aij(ronly,2), dwi(ronly,2) - rs(ronly,1) .* Kij(ronly)));
            end
            
            % s is constrained only
            if ~isempty(sonly)
                rs(sonly,2) = min(this.C - aij(sonly,2),max(-aij(sonly,2), rs(sonly,2)));
                rs(sonly,1) = min(this.C-aij(sonly,1), max(-aij(sonly,1), dwi(sonly,1) - rs(sonly,2) .* Kij(sonly)));
            end
            
            % Case of both vars constrained
            % Update s2 and r2 for the cases of the other one is constrained
            if ~isempty(both)
                Kijb = Kij(both);
                r1 = min(this.C - aij(both,1),max(- aij(both,1), rs(both,1)));
                s1 = min(this.C-aij(both,2), max(-aij(both,2), dwi(both,1) - r1 .* Kijb));
                s2 = min(this.C - aij(both,2),max(- aij(both,2), rs(both,2)));
                r2 = min(this.C-aij(both,1), max(-aij(both,1), dwi(both,2) - s2 .* Kijb));
                
                % Determine update via get higher gain
                % Gain if r is constrained and s accordingly updated
                g1 = r1 .* (dwi(both,1) -.5*r1) + s1 .*(dwi(both,2) -.5 * s1) + sgn(both) .* r1 .* s1.* Kijb;
                % Gain if s is constrained and r accordingly updated
                g2 = r2 .* (dwi(both,1) -.5*r2) + s2 .*(dwi(both,2) -.5 * s2) + sgn(both) .* r2 .* s2 .* Kijb;

                %% See which gain is higher.
                cmp = g1 > g2;
                rs(both(cmp),1) = r1(cmp);
                rs(both(cmp),2) = s1(cmp);
                rs(both(~cmp),1) = r2(~cmp);
                rs(both(~cmp),2) = s2(~cmp);
            end
            
            %% Compute gain
            g = sum(rs .* (dwi - .5*rs),2) + sgn.*rs(:,1).*rs(:,2).*Kij;
            
            % Extract r,s updates
            [dummy, idx] = max(g);
            r = rs(idx,1);
            s = rs(idx,2);
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
            svr.Eps = .1;
            svr.Lambda = 1/20;%1/20; % i.e. C=10 as in ScalarEpsSVR
            svr.Vis = 0;
            
            %kernel = kernels.PolyKernel(7);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(.8);
            svr.K = kernel.evaluate(x,x);

            [ai, svidx] = svr.computeKernelCoefficients(fx,[]);
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
            plot(sv,fx(svidx),'.r','MarkerSize',20);
            plot(x(skipped),fx(skipped),'xr');
            
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