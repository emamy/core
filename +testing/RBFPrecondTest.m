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
        function run
%             seed = 2;
            seed = cputime * 1000;
            rnd = RandStream('mt19937ar','Seed',seed);
            
            
            % Random centers  
            c = 100;
            dim = 300;
            r = [-10 10];
            x = sort(r(1) + rnd.rand(dim,c)*(r(2)-r(1)),2);
            dist = 4e9;

            % Example 2.2 from S08
%             dist = 5e4;
%             r = [0 1];
%             dim = 2;
%             x = [(0:5)/5; (0:5).^2/25];
%             c = size(x,2);
            
            k = kernels.GaussKernel;
            gam = k.setGammaForDistance(dist,eps);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = k;
            kexp.Centers.xi = x;
            
            K = kexp.getKernelMatrix;
            ki = general.interpolation.KernelInterpol;
            ki.init(data.MemoryKernelMatrix(K));
            kexp.Ma = ki.interpolate(f(x))';
            
            pre_kexp = kexp.clone;
            P = getPreconditioner(x, k);
            ki.init(data.MemoryKernelMatrix(P*K));
            pre_kexp.Ma = ki.interpolate(f(x)*P')';
            
            t = PrintTable;
            t.HasHeader = true;
            t.Caption = sprintf('Comparison results for preconditioning test, dim=%d, #c=%d, gamma=%f, rnd-seed=%f',dim,c,gam,seed);
            t.addRow('Results','no preconditioning','with preconditioning','difference');
            cK = cond(K);
            cPK = cond(P*K);
            t.addRow('Condition numbers',cK,cPK,cK-cPK,{'%s','%e','%e','%e'});
            t.addRow('Ma_norms',sum(kexp.Ma_norms),sum(pre_kexp.Ma_norms),...
                sum(kexp.Ma_norms)-sum(pre_kexp.Ma_norms),{'%s','%e','%e','%e'});
            
            % Error on centers
            valfx = f(x);
            avalfx = kexp.evaluate(x);
            pavalfx = pre_kexp.evaluate(x);
            caerr = sum(abs(valfx-avalfx));
            cpaerr = sum(abs(valfx-pavalfx));
            
            % Error on random set
            valx = r(1) + rnd.rand(dim,10000)*(r(2)-r(1));
            valfx = f(valx);
            avalfx = kexp.evaluate(valx);
            pavalfx = pre_kexp.evaluate(valx);
            aerr = sum((valfx-avalfx).^2);
            paerr = sum((valfx-pavalfx).^2);
            
            t.addRow('Error on centers',caerr,cpaerr,caerr-cpaerr,{'%s','%e','%e','%e'});
            t.addRow('Error on random validation set',aerr,paerr,aerr-paerr,{'%s','%e','%e','%e'});
            
            t.print;
            
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
            
            function [P, PM] = getPreconditioner(x, k)
                N = size(x,2);
                M = [];
                I_N = eye(N);
                pivi = 1;
                cols = 1;
                L = I_N; P = I_N;
                tk = zeros(1,N);
                mi = general.MonomialIterator(dim);
                
                while pivi <= N
                    % Compute next monomial
                    
                    % Stragegy 1: Random
%                     deg = ceil(randn(1)*N);
%                     alpha = mi.getRandMonomial(deg);
                    
                    % Stragegy 2: Ordered list
                    if pivi == 1
                        alpha = mi.getNullMonomial;
                        deg = 0;
                    else
                        [alpha, deg] = mi.nextMonomial;
                    end         
                    
                    % Compute new moment matrix column
                    newMcol = prod(x .^ repmat(alpha,1,N),1)';
                    M = [M newMcol];%#ok
                    
                    % apply previous changes to new column
                    newMcol = L*P*newMcol;
                    
                    % Select pivot element candidate indices in current column
                    sel = pivi:N;

                    % Strategy one: Use maximum pivoting
%                     [v, maxidx] = max(abs(newMcol(sel)));
%                     permidx = sel(maxidx);
                    
                    % Strategy two: Only find first nonzero-row and use it
                    permidx = sel(find(abs(newMcol(sel)) > sqrt(eps),1));
                    if ~isempty(permidx)
                        v = abs(newMcol(permidx));
                    else
                        v = 0;
                    end
                    
                    % step one column ahead if current column is already annihilated
                    if v < sqrt(eps)
                        %v
                        cols = cols+1;
                        continue;
                    end
                    
                    % get new permutation matrix according to pivot
                    Pn = getPermMat(pivi, permidx, N);
                    % swap columns
                    newMcol = Pn*newMcol;
                    % compose L_k
                    Ln = I_N;
                    l_k = -newMcol(pivi+1:N)/newMcol(pivi);
                    Ln(pivi+1:N,pivi) = l_k;
                    
                    % keep record of the column indices at which the next linear independent
                    % monomial was added
                    tk(pivi) = deg;
                    
                    L = Ln*Pn*L*Pn;
                    P = Pn*P;
                    
                    pivi = pivi+1;
                    cols = cols+1;
                end
                
                %U = L*P*M;
                D = diag(k.Gamma.^tk);
                PM = P; % return accum. permutation matrix
                P = D * L* P; % compute preconditioning matrix
                
                function P = getPermMat(i,j,n)
                    P = eye(n);
                    P(i,i) = 0; P(j,j) = 0;
                    P(j,i) = 1; P(i,j) = 1;
                end
            end
        end
    end
end