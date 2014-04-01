classdef MatUtils
    %MatUtils: Matrix utility functions
    %
    % @new{0,1,dw,2010-06-02} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing    
    
    methods(Static)
        function [A, idxmat] = laplacemat(h, d1, d2)
            %Computes a 2D diffusion sparse matrix with zero neuman
            %boundary conditions.
            %
            % Arguments:
            % h - discretization spatial stepwidth
            % d1 - region dimension 1
            % d2 - region dimension 2
            %
            % @author Daniel Wirtz @date 11.03.2010
            
            m = d1*d2;
            idxmat = zeros(d1,d2);
            idxmat(:) = 1:m;
            
            inner = find(conv2(ones(size(idxmat)),[0 1 0;1 0 1; 0 1 0],'same')==4);
            
            i = []; j = []; s =[];
            %% Inner points
            % N E S W
            createStencil([-1 d1 1 -d1],[-4 1 1 1 1], inner);
            
            %% Points that are on a boundaries
            % Left boundary neumann
            createStencil([-1 d1 1] ,[-3 1 1 1],(2:d1-1)');
            % Right boundary neumann
            createStencil([-1 -d1 1] ,[-3 1 1 1],(m-d1+2:m-1)');
            % Top boundary neumann
            createStencil([-d1 1 d1] ,[-3 1 1 1],(d1+1:d1:m-2*d1+1)');
            % Bottom boundary neumann
            createStencil([-d1 -1 d1] ,[-3 1 1 1],(2*d1:d1:m-d1)');
            
            % Edge points
            % Top left
            createStencil([1 d1] ,[-2 1 1],1);
            % Top right
            createStencil([1 -d1] ,[-2 1 1],m-d1+1);
            % bottom left
            createStencil([-1 d1] ,[-2 1 1],d1);
            % bottom right
            createStencil([-1 -d1] ,[-2 1 1],m);
            %% Compile matrix
            A = (1/h^2)*sparse(i,j,s,m,m);
            
            % Some comment on how createStencil works - ABOVE
            function createStencil ( stencil , weights , points )
                % Some comment on how createStencil works - BELOW
                %
                % Parameters:
                % stencil - Testdescription
                %
                % weights - bla bla
                
                i = [i; idxmat(points)];
                j = [j; idxmat(points)];
                s = [s; weights(1)* ones(size(points))];
                idx=2;
                for offset = stencil
                    i = [i; idxmat(points)];
                    j = [j; idxmat(points+offset)];
                    s = [s; weights(idx)*ones(size(points))];
                    idx=idx+1;
                end
            end
        end
        
        function A = laplacemat3D(h, g)
            %Computes a 3D laplacian sparse matrix
            %
            % Arguments:
            % h: discretization spatial stepwidth
            % g: the general.geometry.RectGrid3D instance describing the
            % geometry.
            %
            % @author Daniel Wirtz @date 2011-10-06
            
            i = []; j = []; s =[];
            m = g.Dims(1);
            k = g.Dims(2);
            L = g.IndexMatrix;
            
            %% Inner points
            createStencil([-1 m 1 -m -m*k m*k],[6 -1 -1 -1 -1 -1 -1], g.Inner);
            
            %% Points that are on inner boundaries
            % Front boundary neumann
            createStencil([-1 m 1 -m m*k] ,[5 -1 -1 -1 -1 -1],g.Sides.F);
            % Back boundary neumann
            createStencil([-1 m 1 -m -m*k] ,[5 -1 -1 -1 -1 -1],g.Sides.Ba);
            % Left boundary neumann
            createStencil([-1 1 m m*k -m*k] ,[5 -1 -1 -1 -1 -1],g.Sides.L);
            % Right boundary neumann
            createStencil([-1 1 -m m*k -m*k] ,[5 -1 -1 -1 -1 -1],g.Sides.R);
            % Top boundary neumann
            createStencil([1 m -m m*k -m*k] ,[5 -1 -1 -1 -1 -1],g.Sides.T);
            % Bottom boundary neumann
            createStencil([-1 m -m m*k -m*k] ,[5 -1 -1 -1 -1 -1],g.Sides.Bo);

            % Corner points
            %% Front side corners
            % Front top
            createStencil([1 m -m m*k] ,[4 -1 -1 -1 -1],g.Corners.FT);
            % Front bottom
            createStencil([-1 m -m m*k] ,[4 -1 -1 -1 -1],g.Corners.FBo);
            % Front left
            createStencil([1 -1 m m*k] ,[4 -1 -1 -1 -1],g.Corners.FL);
            % Front right
            createStencil([1 -1 -m m*k] ,[4 -1 -1 -1 -1],g.Corners.FR);
            
            %% Back side corners
            % Back top
            createStencil([1 m -m -m*k] ,[4 -1 -1 -1 -1],g.Corners.BaT);
            % Back bottom
            createStencil([-1 m -m -m*k] ,[4 -1 -1 -1 -1],g.Corners.BaBo);
            % Back left
            createStencil([1 -1 m -m*k] ,[4 -1 -1 -1 -1],g.Corners.BaL);
            % Back right
            createStencil([1 -1 -m -m*k] ,[4 -1 -1 -1 -1],g.Corners.BaR);
            
            %% front-back connecting corners
            % top left
            createStencil([1 m m*k -m*k] ,[4 -1 -1 -1 -1],g.Corners.TL);
            % top right
            createStencil([1 -m m*k -m*k] ,[4 -1 -1 -1 -1],g.Corners.TR);
            % bottom left
            createStencil([-1 m m*k -m*k] ,[4 -1 -1 -1 -1],g.Corners.BoL);
            % bottom right
            createStencil([-1 -m m*k -m*k] ,[4 -1 -1 -1 -1],g.Corners.BoR);

            %% Edge points
            % Front Top left
            createStencil([1 m m*k] ,[3 -1 -1 -1],g.Edges.FTL);
            % Front Top right
            createStencil([1 -m m*k] ,[3 -1 -1 -1],g.Edges.FTR);
            % Front bottom left
            createStencil([-1 m m*k] ,[3 -1 -1 -1],g.Edges.FBoL);
            % Front bottom right
            createStencil([-1 -m m*k] ,[3 -1 -1 -1],g.Edges.FBoR);
            % Rear Top left
            createStencil([1 m -m*k] ,[3 -1 -1 -1],g.Edges.BaTL);
            % Rear Top right
            createStencil([1 -m -m*k] ,[3 -1 -1 -1],g.Edges.BaTR);
            % Rear bottom left
            createStencil([-1 m -m*k] ,[3 -1 -1 -1],g.Edges.BaBoL);
            % Rear bottom right
            createStencil([-1 -m -m*k] ,[3 -1 -1 -1],g.Edges.BaBoR);
            
            %% Compile matrix
            A = -(1/h^2)*sparse(i,j,s,g.Points,g.Points);
            
            % Some comment on how createStencil works - ABOVE
            function createStencil ( stencil , weights , points )
                % Some comment on how createStencil works - BELOW
                %
                % Parameters:
                % stencil - Testdescription
                %
                % weights - bla bla
                
                i = [i; L(points)];
                j = [j; L(points)];
                s = [s; weights(1)* ones(size(points))];
                cnt=2;
                for offset = stencil
                    i = [i; L(points)]; %#ok
                    j = [j; L(points+offset)];%#ok
                    s = [s; weights(cnt)*ones(size(points))];%#ok
                    cnt=cnt+1;
                end
            end
        end
        
        function A = generalizedLaplacemat3D(h, g, sigma)
            % Computes a 3D generalized laplacian sparse matrix
            % `\nabla\cdot\sigma\nabla`
            %
            % Arguments:
            % h: discretization spatial stepwidth in all three directions
            % g: the general.geometry.RectGrid3D instance describing the
            % geometry.
            % sigma: the conductivity tensor, assumed to be diagonal
            %
            % @author Timm Strecker @date 2014-03-20
            
            s1 = sigma(1,1)/h(1)^2; s2 = sigma(2,2)/h(2)^2; s3 = sigma(3,3)/h(3)^2;
            %h1 = h(1)^2; h2 = h(2)^2; h3 = h(3)^2;
            i = []; j = []; s =[];
            m = g.Dims(1);
            k = g.Dims(2);
            
            %% Inner points
            createStencil([-1 m 1 -m -m*k m*k],[2*(s1+s2+s3) -s1 -s2 -s1 -s2 -s3 -s3], g.Inner);
            
            %% Points that are on inner boundaries
            % Front boundary neumann
            createStencil([-1 m 1 -m m*k] ,[2*s1+2*s2+s3 -s1 -s2 -s1 -s2 -s3],g.Sides.F);
            % Back boundary neumann
            createStencil([-1 m 1 -m -m*k] ,[2*s1+2*s2+s3 -s1 -s2 -s1 -s2 -s3],g.Sides.Ba);
            % Left boundary neumann
            createStencil([-1 1 m m*k -m*k] ,[2*s1+s2+2*s3 -s1 -s1 -s2 -s3 -s3],g.Sides.L);
            % Right boundary neumann
            createStencil([-1 1 -m m*k -m*k] ,[2*s1+s2+2*s3 -s1 -s1 -s2 -s3 -s3],g.Sides.R);
            % Top boundary neumann
            createStencil([1 m -m m*k -m*k] ,[s1+2*s2+2*s3 -s1 -s2 -s2 -s3 -s3],g.Sides.T);
            % Bottom boundary neumann
            createStencil([-1 m -m m*k -m*k] ,[s1+2*s2+2*s3 -s1 -s2 -s2 -s3 -s3],g.Sides.Bo);

            % Corner points
            %% Front side corners
            % Front top
            createStencil([1 m -m m*k] ,[s1+2*s2+s3 -s1 -s2 -s2 -s3],g.Corners.FT);
            % Front bottom
            createStencil([-1 m -m m*k] ,[s1+2*s2+s3 -s1 -s2 -s2 -s3],g.Corners.FBo);
            % Front left
            createStencil([1 -1 m m*k] ,[2*s1+s2+s3 -s1 -s1 -s2 -s3],g.Corners.FL);
            % Front right
            createStencil([1 -1 -m m*k] ,[2*s1+s2+s3 -s1 -s1 -s2 -s3],g.Corners.FR);
            
            %% Back side corners
            % Back top
            createStencil([1 m -m -m*k] ,[s1+2*s2+s3 -s1 -s2 -s2 -s3],g.Corners.BaT);
            % Back bottom
            createStencil([-1 m -m -m*k] ,[s1+2*s2+s3 -s1 -s2 -s2 -s3],g.Corners.BaBo);
            % Back left
            createStencil([1 -1 m -m*k] ,[2*s1+s2+s3 -s1 -s1 -s2 -s3],g.Corners.BaL);
            % Back right
            createStencil([1 -1 -m -m*k] ,[2*s1+s2+s3 -s1 -s1 -s2 -s3],g.Corners.BaR);
            
            %% front-back connecting corners
            % top left
            createStencil([1 m m*k -m*k] ,[s1+s2+2*s3 -s1 -s2 -s3 -s3],g.Corners.TL);
            % top right
            createStencil([1 -m m*k -m*k] ,[s1+s2+2*s3 -s1 -s2 -s3 -s3],g.Corners.TR);
            % bottom left
            createStencil([-1 m m*k -m*k] ,[s1+s2+2*s3 -s1 -s2 -s3 -s3],g.Corners.BoL);
            % bottom right
            createStencil([-1 -m m*k -m*k] ,[s1+s2+2*s3 -s1 -s2 -s3 -s3],g.Corners.BoR);

            %% Edge points
            % Front Top left
            createStencil([1 m m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.FTL);
            % Front Top right
            createStencil([1 -m m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.FTR);
            % Front bottom left
            createStencil([-1 m m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.FBoL);
            % Front bottom right
            createStencil([-1 -m m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.FBoR);
            % Rear Top left
            createStencil([1 m -m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.BaTL);
            % Rear Top right
            createStencil([1 -m -m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.BaTR);
            % Rear bottom left
            createStencil([-1 m -m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.BaBoL);
            % Rear bottom right
            createStencil([-1 -m -m*k] ,[s1+s2+s3 -s1 -s2 -s3],g.Edges.BaBoR);
            
            %% Compile matrix
            A = -sparse(i,j,s,g.Points,g.Points);
            
            % Some comment on how createStencil works - ABOVE
            function createStencil ( stencil , weights , points )
                % Some comment on how createStencil works - BELOW
                %
                % Parameters:
                % stencil - Testdescription
                %
                % weights - bla bla
                
                i = [i; points];
                j = [j; points];
                s = [s; weights(1)* ones(size(points))];
                cnt=2;
                for offset = stencil
                    i = [i; points]; %#ok
                    j = [j; points+offset];%#ok
                    s = [s; weights(cnt)*ones(size(points))];%#ok
                    cnt=cnt+1;
                end
            end
        end
        
        function Kinv = getExtendedInverse(Kinv, v, mode)
            % For a matrix K with known inverse Kinv the matrix is extended by the vector
            % @docupdate
            if nargin < 3
                % Experiments show that mode 3 is best for gaussian kernel matrices
                mode = 3;
            end
            oldn = size(v,1)-1;
            vold = v(1:oldn);
            e = zeros(size(v));
            e(end) = 1;
            v(end) = (v(end)-1)/2;
            U = [v e]; 
            V = [e v]';
            tmp = [Kinv zeros(oldn,1); zeros(1,oldn) 1];
            
            kh = vold'*Kinv*vold;
            det = (1 + v(end) - kh);
            if mode == 1
                hlp = [-vold*vold' vold; vold' 2*v(end)-vold'*Kinv*vold]/det;
                Kinv = tmp*(eye(oldn+1)-hlp*tmp);
                return;
            end
            
            b = (1/det) * [1+v(end) -1; -kh-v(end)^2 1+v(end)];
            if mode == 2
                Kinv = tmp - (tmp*U)*b*(V*tmp);
                return;
            end
            
            if mode == 3
                a = [Kinv*vold zeros(oldn,1); v(end) 1];
                c = [zeros(1,oldn) 1; vold'*Kinv v(end)];
                Kinv = tmp - a*b*c;
            end
        end
        
        function Asparse = toSparse(A, rowidx, n)
            % Converts a full matrix A to a sparse matrix, where the entries of A are split up
            % to the row indices specified by rowidx.
            %
            % Parameters:
            % A: The matrix to "sparsify" @type matrix<double>
            % rowidx: The row indices of A rows inside the sparse matrix @type rowvec<integer>
            % n: The row size of the target sparse matrix.
            o = size(A,2);
            i = repmat(rowidx,o,1);
            j = reshape(repmat(1:o,length(rowidx),1),[],1);
            Asparse = sparse(i,j,A(:),n,o);
        end
    end
    
    methods(Static)
        function test_ExtendedInverseDirect
            % Tests direct inversion extension for a matrix of size s
            % for a fixed number of experiments each.
            %
            % Repeats all experiments for different kernel widths.
            
            % Matrix size
            s = 100;
            
            % Experiments
            n = 30;
            
            % Kernel widths
            w = 20;
            
            g = kernels.GaussKernel();
            gammas = linspace(sqrt(s)/10,3*sqrt(s),w);
            err = zeros(3,n*w);
            for gam=1:w
                g.Gamma = gammas(gam);
                %fprintf('Using gamma = %f...\n',g.Gamma);
                for exp=1:n
                    X = rand(s,s);
                    A = g.evaluate(X, X);
                    y = rand(s,1);
                    v = g.evaluate(X, y);
                    v(end+1) = g.evaluate(y,y);%#ok
                    Bi = inv([A v(1:end-1); v']);
                    i = (gam-1)*n + exp;
                    Ai = inv(A);
                    err(1,i) = errfun(Bi,MatUtils.getExtendedInverse(Ai, v, 1));
                    err(2,i) = errfun(Bi,MatUtils.getExtendedInverse(Ai, v, 2));
                    err(3,i) = errfun(Bi,MatUtils.getExtendedInverse(Ai, v, 3));
                end
            end
            e1best = sum(err(1,:) < min(err(2:3,:),[],1));
            e2best = sum(err(2,:) < min(err([1 3],:),[],1));
            e3best = sum(err(3,:) < min(err(1:2,:),[],1));
            semilogy(1:n*w,err);
            title(sprintf('Plot for %d gamma values with matrix size %d, %d experiments each.\nBest inverse (norm-difference wise): 1:%d, 2:%d, 3:%d',w,s,n,e1best,e2best,e3best));
            axis tight;
            
            function e = errfun(A,B)
                %e = norm(A-B);
                e = max(abs(A(:)-B(:)));
            end
        end
        
        function test_ExtendedInverseSequential
            % 
            
            % Matrix size
            s = 200;
            
            % Kernel widths
            w = 6;
            
            g = kernels.GaussKernel();
            gammas = linspace(sqrt(s)/10,4*sqrt(s),w);
            err = zeros(3,w*(s-1));
            T = zeros(4,w*(s-1));
            X = rand(s,s);
            for gam=1:w
                g.Gamma = gammas(gam);
                %fprintf('Using gamma = %f...\n',g.Gamma);
                Am = g.evaluate(X,X);
%                 fx = rand(s,1);
                Ai1 = 1/Am(1,1);
                Ai2 = 1/Am(1,1);
                Ai3 = 1/Am(1,1);
                for part=1:s-1
                    i = (gam-1)*(s-1) + part;
                    %A = Am(1:part,1:part);
                    v = Am(1:part+1,part+1);
                    t = tic;
                    Ai1 = MatUtils.getExtendedInverse(Ai1, v, 1);
                    T(1,i) = toc(t);
                    t = tic;
                    Ai2 = MatUtils.getExtendedInverse(Ai2, v, 2);
                    T(2,i) = toc(t);
                    t = tic;
                    Ai3 = MatUtils.getExtendedInverse(Ai3, v, 3);
                    T(3,i) = toc(t);
                    t = tic;
                    Bi = inv(Am(1:part+1,1:part+1));
                    T(4,i) = toc(t);
                    
                    err(1,i) = errfun(Bi,Ai1);
                    err(2,i) = errfun(Bi,Ai2);
                    err(3,i) = errfun(Bi,Ai3);
                end
            end
            e1best = sum(err(1,:) < min(err(2:3,:),[],1));
            e2best = sum(err(2,:) < min(err([1 3],:),[],1));
            e3best = sum(err(3,:) < min(err(1:2,:),[],1));
            semilogy(1:w*(s-1),err);
            T = sum(T,2);
            title(sprintf(['Plot for %d gamma values with matrix size %d.\n'...
                'Count for best mode (norm-difference wise): 1:%d, 2:%d, 3:%d\n'...
                'Times for modes: 1: %1.3fs, 2: %1.3fs, 3: %1.3fs, direct: %1.3fs'],w,s,e1best,e2best,e3best,...
                T(1),T(2),T(3),T(4)));
            axis tight;
            
            
            function e = errfun(A,B)
                %e = norm(A-B);
                e = max(abs(A(:)-B(:)));
            end
        end
    end
    
end

