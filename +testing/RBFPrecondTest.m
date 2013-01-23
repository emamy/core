classdef RBFPrecondTest
% RBFPrecondTest: Test for applicability of the preconditioning method described in [S08]:
% R. Schaback, Limit problems for interpolation by analytic radial basis functions,
% J. Comp. Appl. Math., 2008, 212, 127 - 149
%
% @author Daniel Wirtz @date 2012-01-18
%
% @new{0,6,dw,2012-01-18} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        function [condn, err, epsi, man] = run(dist, conf)
            
            if nargin < 2
                conf.seed = cputime * 1000;
                if nargin < 1
                    dist = 4e9;
                end
            end
            rnd = RandStream('mt19937ar','Seed',conf.seed);
            if ~isfield(conf,'estr')
                conf.estr = 'L_{inf}-L_1';
            	conf.efun = @(x)max(sum(abs(x),1));
            end
            if ~isfield(conf,'rel')
                conf.rel = true;
            end
            if ~isfield(conf,'c')
                conf.c = 100;
            end
            if ~isfield(conf,'dim')
                conf.dim = 300;
            end
            
            r = 2*[-10 10];
            x = sort(r(1) + rnd.rand(conf.dim,conf.c)*(r(2)-r(1)),2);

            % Example 2.2 from S08
%             dist = 5e4;
%             r = [0 1];
%             dim = 2;
%             x = [(0:5)/5; (0:5).^2/25];
%             c = size(x,2);
            
            k = kernels.GaussKernel;
            gam = k.setGammaForDistance(dist,eps);
            epsi = k.epsilon;
            kexp = kernels.KernelExpansion;
            kexp.Kernel = k;
            kexp.Centers.xi = x;
            
            ki = general.interpolation.KernelInterpol;
            ki.UsePreconditioning = false;
            ki.init(kexp);
            kexp.Ma = ki.interpolate(f(x))';
            pre_kexp = kexp.clone;
            
            ki.UsePreconditioning = true;
            ki.init(kexp);
            P = ki.getPreconditioner(k, x);
            pre_kexp.Ma = ki.interpolate(f(x)*P')';
            
            K = kexp.getKernelMatrix;
            condn(1) = cond(K);
            condn(2) = cond(P*K);
            
            % Error on centers
            valfx = f(x);
            avalfx = kexp.evaluate(x);
            pavalfx = pre_kexp.evaluate(x);
            err(1) = conf.efun(valfx-avalfx);
            err(2) = conf.efun(valfx-pavalfx);
            
            if conf.rel
                err = err ./ conf.efun(valfx);
            end
            
            man(1) = sum(kexp.Ma_norms);
            man(2) = sum(pre_kexp.Ma_norms);
            if nargout == 0
                t = PrintTable;
                t.HasHeader = true;
                t.Caption = sprintf('Comparison results for preconditioning test, dim=%d, #c=%d, gamma=%f, rnd-seed=%f',conf.dim,conf.c,gam,conf.seed);
                t.addRow('Results','no preconditioning','with preconditioning','difference');
                t.addRow('Condition numbers',condn(1),condn(2),condn(1)-condn(2),{'%s','%e','%e','%e'});
                t.addRow('Ma_norms',man(1),man(2),...
                    man(1)-man(2),{'%s','%e','%e','%e'});
                
                % Error on random set
                valx = r(1) + rnd.rand(conf.dim,10000)*(r(2)-r(1));
                valfx = f(valx);
                avalfx = kexp.evaluate(valx);
                pavalfx = pre_kexp.evaluate(valx);
                aerr = conf.efun(valfx-avalfx);
                paerr = conf.efun(valfx-pavalfx);
                
                if conf.rel
                    aerr = aerr / conf.efun(valfx);
                    paerr = paerr / conf.efun(valfx);
                end
                re = '';
                if conf.rel
                    re = 'Rel. ';
                end
                t.addRow(sprintf('%s%s-Error on centers',re,conf.estr),err(1),err(2),err(1)-err(2),{'%s','%e','%e','%e'});
                t.addRow(sprintf('%s%s-Error on random validation set',re,conf.estr),aerr,paerr,aerr-paerr,{'%s','%e','%e','%e'});
                t.print;
            end
            
%             if dim == 1
%                 xv = repmat(linspace(r(1),r(2),200),dim,1);
%                 plot(xv(1,:),f(xv),'r',xv(1,:),kexp.evaluate(xv)','b',...
%                     xv(1,:),pre_kexp.evaluate(xv)','g');
%             else
%                 atd=data.ApproxTrainData(x,[],[]);
%                 atd.fxi = f(x);
%                 FunVis2D(kexp,atd,[],@f);
%                 FunVis2D(pre_kexp,atd,[],@f);
%             end
            
            function fx = f(x,~,~)
                fx = sum(sin(pi*x/5),1);
            end
        end
        
        function [C, E, M, ep] = runForDists(minr, maxr, c)
            num = 25;
            
            % Random seed
            conf.seed = 2;
            % relative error
            conf.rel = true;
            % centers
            if nargin < 3
                conf.c = 200;
            else
                if c < 2
                    error('Minimum two centers required.');
                end
                conf.c = c;
            end
            % dims
            conf.dim = 1000;
            
            % L1 error
            conf.estr = 'L_{inf}-L_1';
            conf.efun = @(x)max(sum(abs(x),1));
            
            % L2 error
%             conf.estr = 'L_2-L_2';
%             conf.efun = @(x)norm(sqrt(sum(x.^2,1)));
            
            % Linf error
%             conf.estr = 'L_inf';
%             conf.efun = @(x)max(max(abs(x),[],1));
            
            if nargin == 0
                minr = 1e-2;
                maxr = 1e12;
            end
            dist = logspace(log10(minr),log10(maxr),num);
            C = zeros(2,num);
            E = C;
            M = C;
            ep = zeros(1,num);
            fprintf('Starting %d runs for %d centers with dist range from %f to %f, seed=%e...\n',num,conf.c,minr,maxr,conf.seed)
            for i = 1:num
                [C(:,i), E(:,i), ep(i), M(:,i)] = testing.RBFPrecondTest.run(dist(i), conf);
            end
            if nargout == 0
                figure;
                subplot(1,2,1);
                loglog(ep,C);
                legend('no preconditioning','with preconditioning');
                xlabel('epsilon (kernel dilation param)'); ylabel('condition number');
                title(sprintf('Condition number comparison (minr=%e, maxr=%e)',minr,maxr));
                axis tight;
                subplot(1,2,2);
                %semilogx(ep,E);
                loglog(ep,E);
                legend('no preconditioning','with preconditioning');
                re = '';
                if conf.rel
                    re = 'Rel. ';
                    hold on;
                    loglog(ep,1,'k');
                    hold off;
                end
                title(sprintf('%s%s-error on %d centers, dim=%d, ',re,conf.estr,conf.c,conf.dim));
                xlabel('epsilon (kernel dilation param)'); ylabel('L2 error over centers');
                axis tight;
            end
        end
        
        function ep = runForCenters(cmax)
            if nargin == 0
                cmax = 250;
            end
            
            minr = 1e-2;
            maxr = 1e12;
            all = 2:5:cmax;
            ep = zeros(1,cmax);
            for cnum = all
                [C, ~, ~, epsi] = testing.RBFPrecondTest.runForDists(minr, maxr, cnum);
                % The first where preconditioning is better
                pidx = find(C(2,:) < C(1,:),1);
                if isempty(pidx)
                    pidx = size(C,2);
                end
                ep(cnum) = epsi(pidx);
            end
            if nargout == 0
                semilogy(all,ep(all));
                xlabel('Number of centers'); ylabel('epsilon at which preconditioning gives better conditioning');
            end
        end
    end
end