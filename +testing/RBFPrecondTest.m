classdef RBFPrecondTest
% RBFPrecondTest: 
%
%
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
    
    properties
    end
    
    methods(Static)
        function run
            
            % Random centers
%             c = 10;
%             dim = 1;
%             r = [-10 10];
%             x = sort(r(1) + rnd.rand(dim,c)*(r(2)-r(1)),2);
            % Example 2.2 from S08
            c = 6;
            dim = 2;
            x = [(0:5)/5; (0:5).^2/25];
            
            k = kernels.GaussKernel(50000);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = k;
            kexp.Centers.xi = x;
            
            K = kexp.getKernelMatrix;
%             gam = k.setGammaForDistance(20,eps);
%             while true
%                 K = kexp.getKernelMatrix;
%                 if rank(K) < c
%                     k.Gamma = gam*2/3;
%                     K = kexp.getKernelMatrix;
%                     break;
%                 end
%                 gam = gam * 1.5;
%                 k.Gamma = gam;
%             end
%             imagesc(K);

            ki = general.interpolation.KernelInterpol;
            ki.init(data.MemoryKernelMatrix(K));
            kexp.Ma = ki.interpolate(f(x))';
            
            pre_kexp = kexp.clone;
            %P = eye(c); P(5,3) = 0.0001;
            P = getPreconditioner(x, k);
            ki.init(data.MemoryKernelMatrix(P*K));
            pre_kexp.Ma = ki.interpolate(f(x)*P')';
            
            if dim == 1
                xv = repmat(linspace(r(1),r(2),200),dim,1);
                plot(xv(1,:),f(xv),'r',xv(1,:),kexp.evaluate(xv)','b',...
                    xv(1,:),pre_kexp.evaluate(xv)','g');
            else
                atd=data.ApproxTrainData(x,[],[]);
                atd.fxi = f(x);
                FunVis2D(kexp,atd,[],@f);
                FunVis2D(pre_kexp,atd,[],@f);
            end
            
            function fx = f(x,~,~)
                fx = sum(sin(pi*x/5),1);
            end
            
            function P = getPreconditioner(x, k)
                N = size(x,2); % k_2 = N-1
                M = [];
%                 mon = [];
                I_N = eye(N);
                pivi = 1;
                cols = 1;
                L = I_N; P = I_N;
%                 Li = I_N;
                tk = zeros(1,N);
                mi = general.MonomialIterator(dim);
                
%                 M = zeros(N,N);
%                 M(:,1) = 1;
%                 for i=2:300
%                     % get new raw M column
%                     alpha = repmat(mi.nextMonomial,1,N);
%                     newMcol = prod(x .^ alpha)';
%                     M(:,i) = newMcol;
%                 end
                
                while pivi <= N
                    % Compute next monomial
                    %alpha = mi.getRandMonomial(rnd.randi(N));
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
%                     newMcol = L*P*M(:,cols);
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
                
%                 for i = 1:nd
%                     alpha = repmat(getMonomial(degs(i), dim),1,N);
%                     M(:,i) = sum(x .^ alpha)';
%                 end
%                 [l,u] = lu(M);
%                 l = inv(L);
%                 u = L*P*M;
                U = L*P*M
%                 Li = inv(L);
%                 P*M
                D = diag(k.Gamma.^tk);
                P = D * L* P;
                
                function P = getPermMat(i,j,n)
                    P = eye(n);
                    P(i,i) = 0; P(j,j) = 0;
                    P(j,i) = 1; P(i,j) = 1;
                end
            end
        end
    end
end