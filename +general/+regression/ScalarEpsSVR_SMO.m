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
    end
    
    methods
        function ai = regress(this, fxi)
            
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
                    dW(1:n) = dW(1:n) - r * Ki;
                else
                    dW(n+1:end) = dW(n+1:end) - r * Ki;
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
        
        function ai = regress2(this, fxi)
            
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
            svr.Lambda = eps;%1/20; % i.e. C=10 as in ScalarEpsSVR
            
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