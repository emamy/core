classdef Speed
    % Speed: Collects tests regarding speed of different methods and strategies
    %
    % @author Daniel Wirtz @date 2013-08-30
    %
    % @new{0,7,dw,2013-08-30} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties
    end
    
    methods(Static)
        
        function [times, e] = KernelExpCustomBaseEval(kexp, numpts)
            % Tests the evaluation speed and determines the evaluation
            % error of a kernel expansion and this expansion using the
            % default direct translate base.
            %
            % Must have a custom base set (HasCustomBase = true)
            %
            % Parameters:
            % kexp: A kernel expansion @type kernels.KernelExpansion
            % numpts: The number of random points at which to evaluate @type integer @default
            % 1000
            %
            % Return values:
            % times: A `2\times 1` vector containing the evaluation times using the
            % custom and direct base in the first and second entry, respectively. @type
            % colvec<double>
            % e: The maximum pointwise absolute and relative `L^2`-errors over all runs @type
            % double
            if nargin < 2
                numpts = 1000;
                if nargin < 1
                    kexp = kernels.KernelExpansion;
                    kexp.Centers.xi = rand(400,40);
                    kexp.Ma = rand(30,40);
                    o = general.Orthonormalizer;
                    A = o.orthonormalize(rand(40,40));
                    L = chol(A'*diag(randi(10,1,40))*A);
                    kexp.Base = L;
                end
            end
            if ~kexp.HasCustomBase
                error('No custom base set to compare to.');
            end
            x = rand(size(kexp.Centers.xi,1),numpts);
            dbase = kexp.toTranslateBase;
            runs = 10;
            t = zeros(2,runs);
            e = zeros(2,runs);
            for r = 1:runs
                tic;
                fx1 = kexp.evaluate(x);
                t(1,r) = toc;
                tic;
                fx2 = dbase.evaluate(x);
                t(2,r) = toc;
                er = Norm.L2(fx1-fx2);
                e(1,r) = max(er);
                fxin = Norm.L2(fx1);
                fxin(fxin == 0) = 1;
                e(2,r) = max(er./fxin);
            end
            times = mean(t,2);
            e = max(e,[],2);
        end
        
        function pt = TryCatch(loopsize)
            % TryCatch: Demonstrate how slow try-catch blocks are.
            %
            % Copied from http://www.mathworks.com/matlabcentral/newsreader/view_thread/275243.
            %
            % @author Daniel Wirtz @date 2011-05-18
            %
            % @new{0,4,dw,2011-05-18} Added this function.
            
            t = [];
            
            if nargin < 1
                loopsize = 1000000;
            end
            
            tic
            for ii = 1:loopsize
                A = ii;
            end
            t(end+1) = toc;
            
            tic
            try
                for ii = 1:loopsize
                    A = ii.^2;
                end
            catch ME
            end
            t(end+1) = toc;
            
            tic
            for ii = 1:loopsize
                try
                    A = ii.^2;
                catch ME
                end
            end
            t(end+1) = toc;
            
            tic
            for ii = 1:loopsize
                if ~ii
                    % Never here.
                    %     elseif 0
                    %         % Never here.
                    %     elseif mod(ii,2)
                    %         B = ii;
                else
                    C = ii.^2;
                end
            end
            t(end+1) = toc;
            
            pt = PrintTable;
            pt.HasHeader = true;
            pt.HasRowHeader = true;
            pt.addRow('','Plain for','try(for)','for(try)','for(if)');
            tc = num2cell(t);
            pt.addRow('Time',tc{:},{'%6gs'});
            reltimes = num2cell(t./min(t));
            pt.addRow('Rel factor',reltimes{:},{'%3.3g'});
        end
        
        function pt = BinaryvsMatSave(fill, n , m)
            % BinaryvsMatSaveSpeed: Tests the speeds and storage size required for
            % storing double matrices.
            %
            % Parameters:
            % fill: The percentage in [0,1] of which the matrix should consist of
            % non-zero values @type double @default .5
            % n: Row number @type integer @default 2000
            % m: Row number @type integer @default 1000
            %
            % @author Daniel Wirtz @date 2012-07-09
            %
            % @new{0,6,dw,2012-07-09} Added this function.
            %
            % This class is part of the framework
            % KerMor - Model Order Reduction using Kernels:
            % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
            % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
            % - \c License @ref licensing
            
            if nargin < 3
                m = 1000;
                if nargin < 2
                    n = 2000;
                    if nargin < 1
                        fill = .5;
                    end
                end
            end
            
            a = KerMor.App;
            
            runs = 2;
            pt = PrintTable;
            pt.addRow('','Direct','-v4','-v6','binary');
            pt.HasHeader = true;
            pt.HasRowHeader = true;
            for k = 1:runs
                A = full(sprand(n,m,fill));
                
                f = fullfile(a.DataDirectory,'tmp.mat');
                ti = tic;
                save(f,'A');
                t{1} = toc(ti);
                d = dir(f);
                s{1} = d.bytes;
                delete(f);
                
                f = fullfile(a.DataDirectory,'tmp2.mat');
                ti = tic;
                save(f,'A','-v4');
                t{2} = toc(ti);
                d = dir(f);
                s{2} = d.bytes;
                delete(f);
                
                f = fullfile(a.DataDirectory,'tmp3.mat');
                ti = tic;
                save(f,'A','-v6');
                t{3} = toc(ti);
                d = dir(f);
                s{3} = d.bytes;
                delete(f);
                
                f = fullfile(a.DataDirectory,'tmp4.mat');
                ti = tic;
                fh = fopen(f,'w+');
                fwrite(fh,A,'double');
                fclose(fh);
                t{4} = toc(ti);
                d = dir(f);
                s{4} = d.bytes;
                delete(f);
                pt.addRow('Time',t{:},{'%6gs'});
                pt.addRow('Size',s{:},{@(v)sprintf('%5gk',v/1024)});
            end
        end
        
        function FindVecInMatrix(n,m)
            % Created  for test purposes of finding a vector in a matrix.
            %
            % Parameters:
            % n: Number of matrix rows @type integer @default 500
            % m: Number of matrix columns @type integer @default 6000
            %
            % @author Daniel Wirtz @date 2011-04-01
            %
            % @new{0,5,dw,2011-04-01} Added this function.
            
            if nargin < 2
                m = 12000;
                if nargin < 1
                    n = 5000;
                end
            end
            A = rand(n,m);
            pos = [1 round(m/2) m];
            for kpos = 1:3
                b = A(:,pos(kpos));
                
                time = tic;
                hlp = bsxfun(@eq,A,b);
                i = find(sum(hlp,1) == n,1);
                time1 = toc(time);
                
                time = tic;
                i = find(sum(A == repmat(b,1,m)) == n,1);
                time2 = toc(time);
                
                time = tic;
                ind = strfind(reshape(A,1,[]),b');
                round((ind+n-1)/n);
                time3 = toc(time);
                
                fprintf('Dims: %dx%d, match pos: %d, bsxfun time: %f, repmat time: %f, strfind time: %f\n',...
                    n,m,pos(kpos), time1,time2,time3);
            end
        end
        
        function AffParamMatrix(n)
            % Tests the evaluation speed of the general.AffParamMatrix
            % compared to loop-type evaluations of affine decompositions.
            %
            % Parameters:
            % n: The number of test runs to perform. @type integer
            % @default 5000
            %
            % @author Daniel Wirtz @date 2011-10-25
            %
            % @new{0,5,dw,2011-10-25} Added this function.
            
            if nargin < 1
                n = 5000;
            end
            
            d = 50;
            
            str = {'4*t','4*t + mu(1)','exp(sum(mu))','sin(t)','exp(-t)*mu(2)','t*mu(2)','cos(t)','exp(-mu(2))'};
            theta{1} = @(t,mu)4*t;
            theta{2} = @(t,mu)4*t + mu(1);
            theta{3} = @(t,mu)exp(sum(mu));
            theta{4} = @(t,mu)sin(t);
            theta{5} = @(t,mu)exp(-t)*mu(2);
            theta{6} = @(t,mu)t*mu(2);
            theta{7} = @(t,mu)cos(t);
            theta{8} = @(t,mu)exp(-mu(2));
            nc = length(theta);
            
            am = general.AffParamMatrix;
            
            for i=1:nc
                mat{i} = rand(d,d);%#ok
                am.addMatrix(str{i}, mat{i});
            end
            
            v = rand(d,1);
            params = rand(3,n);
            
            %% Linear
            fprintf('Testing linear case...\n');
            t = tic;
            for i=1:n
                coeff = zeros(1,5);
                for k=1:nc
                    fun = theta{k};
                    coeff(k) = fun(params(1,i),params(2:3,i));
                end
                dummy = lincomb_sequence(mat,coeff)*v;%#ok
            end
            t1 = toc(t);
            
            t = tic;
            for i=1:n
                dummy = am.compose(params(1,i),params(2:3,i))*v;%#ok
            end
            t2 = toc(t);
            
            fprintf('Direct: %gs, AffParamMatrix: %gs\n',t1,t2);
            
            %% Quadratic
            fprintf('Testing quadratic case...\n');
            t = tic;
            for i = 1:nc
                for j = 1:nc
                    mat2{i,j} = mat{i}*mat{j};%#ok
                end
            end
            for i=1:n
                coeff = zeros(1,5);
                for k=1:nc
                    fun = theta{k};
                    coeff(k) = fun(params(1,i),params(2:3,i));
                end
                dummy = lincomb_sequence2(mat2,coeff,coeff)*v;%#ok
            end
            t1 = toc(t);
            
            t = tic;
            am = am*am;
            for i=1:n
                dummy = am.compose(params(1,i),params(2:3,i))*v;%#ok
                %dummy = am(params(1,i),params(2:3,i))*v;%#ok
            end
            t2 = toc(t);
            
            fprintf('Direct: %gs, AffParamMatrix: %gs\n',t1,t2);
            
            function res = lincomb_sequence(seq,sigma)
                %function res = lincomb_sequence(seq,sigma)
                %
                % function performing a linear combination of the elements in the
                % cell array seq with coefficients in sigma result is a
                % vector/matrix the same size as the entries of seq.
                % if sigma = 0, the component is not touched, i.e. the component
                % may also be an empty matrix.
                
                % Bernard Haasdonk 15.5.2007
                
                Q = length(sigma);
                res = seq{1}* sigma(1);
                for q=2:Q
                    if sigma(q)~=0
                        res = res + sigma(q)*seq{q};
                    end;
                end;
            end
            
            function res = lincomb_sequence2(seq,sigma1,sigma2)
                %function res = lincomb_sequence2(seq,sigma1,sigma2)
                %
                % function performing a linear combination of the elements in the
                % 2d cell array seq with coefficients in sigma1 and sigma2 result is a
                % vector/matrix the same size as the entries of seq.
                % size of seq is length(sigma1) x length(sigma2)
                % if some sigma = 0, the component is not touched, i.e. the component
                % may also be an empty matrix.
                
                % Bernard Haasdonk 15.5.2007
                
                Q1 = length(sigma1);
                Q2 = length(sigma2);
                res = zeros(size(seq{1,1}));
                for q1=1:Q1
                    if sigma1(q1)~=0
                        for q2=1:Q2
                            if sigma2(q2)~=0
                                res = res + sigma1(q1)*sigma2(q2)*seq{q1,q2};
                            end;
                        end;
                    end;
                end;
            end
            
        end
        
        function res = GaussMexSpeedTest1Arg(sx,sy,iter)
            % Tests the speed of the c implementations of evaluate for
            % gaussians.
            %
            % One-Argument test (self-evaluation)
            %
            % @note You need to compile/mex the files in the
            % +kernels/@GaussKernel folder for this to work.
            if nargin < 3
                iter = 50;
                if nargin < 2
                    sy = 100;
                    if nargin < 1
                        sx = 5000;
                    end
                end
            end
            k = kernels.GaussKernel(1);            
            x = rand(sx,sy);
            
            fprintf('One argument speed test with sx=%d, sy=%d and %d iterations\n',sx,sy,iter);
            tmex = zeros(1,iter); tmex2 = zeros(1,iter);
            tmexp = zeros(1,iter); tmat = zeros(1,iter);
            pi = ProcessIndicator('Iterating..',iter);
            for i=1:iter
                t = tic;
                Kmex = k.dontuse_evaluate(x);
                %Kmex = k.evaluateIntel(x);
                tmex(i) = toc(t);

                t = tic;
                Kmex2 = k.dontuse_evaluateDirect(x);
                tmex2(i) = toc(t);

                t = tic;
                KmexP = k.evaluateMex(x);
                tmexp(i) = toc(t);

                t = tic;
                K = k.evaluate(x,[]);
                tmat(i) = toc(t);
                pi.step;
            end
            pi.stop;
            fprintf(['1: %1.5fs - Mex straight\n2: %1.5fs - Mex time opt\n'...
                '3: %1.5fs - Mex time opt openmp\n4: %1.5fs - Matlab time\n'...
                'Difference norms: 1-4=%1.5f, 2-4=%1.5f, 3-4=%1.5f\n'],...
                mean(tmex2),mean(tmex),mean(tmexp),mean(tmat),...
                norm(Kmex2-K),norm(Kmex-K),norm(KmexP-K));
            
            res = true;
        end
        
        function res = GaussMexSpeedTest2Arg(sx,sy1,sy2,iter)
            % Tests the speed of the c implementations of evaluate for
            % gaussians.
            %
            % Two-Argument test.
            %
            % @note You need to compile/mex the files in the
            % +kernels/@GaussKernel folder for this to work.
            if nargin < 4
                iter = 40;
                if nargin < 3
                    sy2 = 100;
                    if nargin < 2
                        sy1 = 100;
                        if nargin < 1
                            sx = 500;
                        end
                    end
                end
            end
            
            k = kernels.GaussKernel(1);            
            x = rand(sx,sy1);
            y = rand(sx,sy2);
            
            fprintf('Two argument speed test with sx=%d, sy1=%d, sy2=%d and %d iterations\n',sx,sy1,sy2,iter);
            tmex = zeros(1,iter); tmex2 = zeros(1,iter);
            tmexp = zeros(1,iter); tmat = zeros(1,iter);
            fprintf('Iteration ');
            for i=1:iter
                fprintf('%d ',i);
                
                t = tic;
                Kmex = k.dontuse_evaluate(x,y);
                %Kmex = k.evaluateIntel(x,y);
                tmex(i) = toc(t);

                t = tic;
                Kmex2 = k.dontuse_evaluateDirect(x,y);
                tmex2(i) = toc(t);

                t = tic;
                KmexP = k.evaluateMex(x,y);
                tmexp(i) = toc(t);

                t = tic;
                K = k.evaluate(x,y);
                tmat(i) = toc(t);
            
            end
            fprintf('done!\n');
            fprintf(['1: %1.5fs - Mex straight\n2: %1.5fs - Mex time opt\n'...
                '3: %1.5fs - Mex time opt openmp\n4: %1.5fs - Matlab time\n'...
                'Difference norms: 1-4=%1.5f, 2-4=%1.5f, 3-4=%1.5f\n'],...
                mean(tmex2),mean(tmex),mean(tmexp),mean(tmat),...
                norm(Kmex2-K),norm(Kmex-K),norm(KmexP-K));
            
            res = true;
        end
    end
end