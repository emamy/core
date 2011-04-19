classdef qpInteriorPoint
    %QPINTERIORPOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function copy = clone(this)
            copy = solvers.qp.qpInteriorPoint;
            copy = clone@solvers.qp.BaseQPSolver(this, copy);
        end
    end
    
    methods(Access=protected)
        
        function [p,d,cflag,info] = solve(this,Q,c,lb,ub,A,lbA,ubA,x0)
            error('Not implemented correctly (mess doesnt work!)');
        end
        
        function [ai,b] = quadprog(this, fxi)
            % Own separate derivation of the opt. problem from Smola
            e = 0.05;
            maxIt = 200;
            
            m = size(this.K,1);
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            
            % Problem setup
            Q = T'*this.K*T;
            c = (this.eps*ones(2*m,1) - T'*fxi);
            A = [ones(1,m) -ones(1,m)];
            u = this.C/m;
            
            %% Compute initial guess
            mu = 1;
            I = [Q+diag(ones(size(Q,1),1)) -A'; -A 0];
            initial = I\[-c;0];
            %beta = max(.01*u,min(.9*u,initial(1:end-1)));
            %beta = min(.9*u,initial(1:end-1));
            beta = max(1,initial(1:end-1));
            gamma = initial(end);
            t = max(1,u-beta);
            tau = mu./t;
            lambda = mu./beta;
            gap = Inf;
            cnt = 0;
            
            %% Iterations
            while -log10(gap) < 4 && cnt < maxIt
                % Compose cholesky factorization
                %                 R = chol(Q + diag(tau./t+lambda./beta));
                %
                %                 % Compute helper q vector
                %                 q = (R'\A')';
                %                 % Compose bigger lower and upper triangular matrices
                %                 L = [R' zeros(2*m,1); q 1];
                %                 R = [R q'; zeros(1,2*m) -q*q'];%#ok
                
                % Compose right hand sides
                rhs = [-Q*beta-c-A'*gamma-mu*(1./t + 1./beta);...
                    -A*beta]; %+tau'*(u+beta+tau)./t
                
                % Solve equation system - predictor step
                M = [Q+diag(tau./t+lambda./beta) -A'; -A 0];
                sol = M\rhs;
                %                 sol = R\rhs;
                %                 sol = L\sol;
                dbeta = sol(1:2*m);
                dt = u-beta-t-dbeta;
                dtau = mu./tau.*dt - tau;
                dlambda = mu./lambda.*dbeta - lambda;
                
                % Update rhs side with predicted values
                rhs(1:2*m) = rhs(1:2*m) + dbeta.*dlambda./beta - dt.*dtau./t;
                
                % Solve again - corrector step
                sol = M\rhs;
                %                 sol = R\rhs;
                %                 sol = L\sol;
                % Compute corrected values
                dbeta = sol(1:2*m);
                dt = u-beta-t-dbeta;
                dtau = mu./tau.*dt - tau - dt.*dtau./t;
                dlambda = mu./lambda.*dbeta - lambda - dbeta.*dlambda./beta;
                dgamma = sol(end);
                
                % Compute lamdba factor
                comp = (e-1)*[beta./dbeta t./dt tau./dtau lambda./dlambda];
                nu = min(1,min(comp(comp>=0)));
                
                % Compute new mu
                mu = (tau'*t/(2*m))*((1-nu+e)/(10+nu))^2;
                
                % Update variables
                beta = beta + nu*dbeta;
                t = t + nu*dt;
                tau = tau + nu*dtau;
                lambda = lambda + nu*dlambda;
                gamma = gamma + nu*dgamma;
                
                % Compute cancellation check
                
                gap = tau'*t / abs(.5*beta'*Q*beta+c'*beta+.5*beta'*lambda);
                cnt = cnt+1;
                fprintf('Gap at iteration %d: %f, nu:%2.10f, mu:%f\n',cnt,gap,nu,mu);
            end
            % Return correct coefficients!
            ai = T*beta;
            b = -gamma;
        end
        
        function [ai,b] = quadprog2(this, fxi)
            % Implememtation using no transformation of the `\beta \leq u`
            % to `\beta + t = u, t \geq 0`
            
            e = 0.05;
            maxIt = 100;
            
            m = size(this.K,1);
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            
            % Problem setup
            Q = T'*this.K*T;
            c = (this.eps*ones(2*m,1) - T'*fxi);
            A = [ones(1,m) -ones(1,m)];
            u = this.C/m;
            
            %% Compute initial guess
            mu = .1;
            I = [Q+diag(ones(size(Q,1),1)) -A'; -A 0];
            initial = I\[-c;0];
            %beta = max(.009*u,min(.9*u,initial(1:end-1)));
            %beta = min(.9*u,initial(1:end-1));
            %beta = max(.9*u,initial(1:end-1));
            beta = max(.1,initial(1:end-1));
            gamma = initial(end);
            tau = -mu./(beta-u);
            lambda = mu./beta;
            gap = Inf;
            cnt = 0;
            
            %% Iterations
            while -log10(gap) < 4 && cnt < maxIt
                % Compose cholesky factorization
                %                 D = diag(lambda./beta-tau./(beta-u));
                %                 try
                %                     R = chol(Q + D);
                %                 catch ME
                %                     keyboard;
                %                 end
                %                 % Compute helper q vector
                %                 q = (R'\(-A)')';
                %                 % Compose bigger lower and upper triangular matrices
                %                 L = [R' zeros(2*m,1); q 1];
                %                 R = [R q'; zeros(1,2*m) -q*q'];%#ok
                
                % Compose right hand sides
                rhs = [(-Q*beta -c +A'*gamma +mu*(1./beta - 1./(beta-u))); -A*beta];
                
                % Solve equation system - predictor step
                %                 sol = R\rhs;
                %                 sol = L\sol;
                M = [Q+diag(lambda./beta-tau./(beta-u)) -A'; -A 0];
                sol = M\rhs;
                dbeta = sol(1:2*m);
                dtau = mu./(beta+dbeta-u) - tau;
                dlambda = mu./(beta+dbeta) - lambda;
                
                % Update rhs side with predicted values
                rhs(1:2*m) = rhs(1:2*m) - dbeta.*dlambda./beta + dtau.*dbeta./(beta-u);
                
                % Solve again - corrector step
                %                 sol = R\rhs;
                %                 sol = L\sol;
                sol = M\rhs;
                % Compute corrected values
                dbeta = sol(1:2*m);
                dtau = mu./(beta+dbeta-u) - tau - dtau.*dbeta./(beta-u);
                dlambda = mu./(beta+dbeta) - lambda - dlambda.*dbeta./beta;
                dgamma = sol(end);
                
                % Compute lamdba factor
                comp = (e-1)*[beta./dbeta tau./dtau lambda./dlambda];
                comp = comp(comp>=0);
                if isempty(comp)
                    comp = 1;
                end
                nu = min(1,min(comp));
                
                % Compute new mu
                %mu =
                %(tau'*(beta-u)-lambda'*beta)/(2*m)*((1-nu)/(10+nu))^2;
                mu = (lambda'*beta)/(2*m)*((1-nu+e)/(10+nu))^2;
                
                % Update variables
                beta = beta + nu*dbeta;
                tau = tau + nu*dtau;
                lambda = lambda + nu*dlambda;
                gamma = gamma + nu*dgamma;
                
                % Compute cancellation check
                hlp = tau'*(beta-u)-lambda'*beta;
                %hlp = tau'*(beta-u)-lambda'*beta;
                gap = hlp / (.5*beta'*Q*beta +c'*beta + .5*hlp);
                cnt = cnt+1;
                fprintf('Gap at iteration %d: %f, nu:%f, mu:%f\n',cnt,gap,nu,mu);
            end
            % Return correct coefficients!
            %B = [beta(1:m) beta(m+1:end)];
            %[b,i] = max(abs(B),[],2);
            %ai = b.*sign(B(i));
            
            ai = T*beta;
            b = -gamma;
        end
        
        function [ai,b] = quadprog3(this, fxi)
            % Straight-forward implementation from Smola-Book
            e = 0.05;
            maxIt = 10;
            
            m = size(this.K,1);
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            
            % Problem setup
            Q = T'*this.K*T;
            c = (this.eps*ones(2*m,1) - T'*fxi);
            A = [ones(1,m) -ones(1,m)];
            d = 0;
            u = this.C/m;
            
            %% Compute initial guess
            mu = 1;
            I = [Q+diag(ones(size(Q,1),1)) -A'; -A 0];
            initial = I\[-c;0];
            %alpha = min(1,initial(1:end-1));
            %alpha = max(1,initial(1:end-1));
            alpha = ones(2*m,1)*.0001;
            h = initial(end);
            t = max(1,u-alpha);
            %t = u-alpha;
            s = mu./t;
            z = mu./alpha;
            gap = Inf;
            cnt = 0; nu = Inf;
            
            %% Iterations
            while -log10(gap) < 4 && cnt < maxIt %&& nu > eps
                
                rhop1 = d-A*alpha;
                rhop2 = u-alpha-t;
                rhod = -Q*alpha - c + A'*h - s + z;
                rhokkt1 = mu./alpha - z;
                rhokkt2 = mu./t - s;
                
                % Compose right hand sides
                rhs = [rhod + rhokkt1 - rhokkt2; rhop1];
                
                % Solve equation system - predictor step
                cholflag = true;
                try
                    R = chol(Q + diag(s./t + z./alpha));
                    % Compute helper q vector
                    q = (R'\A')';
                    % Compose bigger lower and upper triangular matrices
                    L = [R' zeros(2*m,1); q 1];
                    R = [R q'; zeros(1,2*m) -q*q'];%#ok
                    sol = R\rhs;
                    sol = L\sol;
                catch ME
                    M = [Q+diag(s./t + z./alpha) -A'; -A 0];
                    sol = M\rhs;
                    cholflag = false;
                end
                dalpha = sol(1:2*m);
                dt = rhop2-dalpha;
                dz = rhokkt1 - (z.*dalpha)./alpha;
                ds = rhokkt2 + (s.*dalpha)./t;
                
                % Update rhs side with predicted values
                rhokkt1 = rhokkt1 - (dalpha.*dz)./alpha;
                rhokkt2 = rhokkt2 - (dt.*ds)./t;
                rhs = [rhod + rhokkt1 - rhokkt2; rhop1];
                
                % Solve again - corrector step
                if cholflag
                    sol = R\rhs;
                    sol = L\sol;
                else
                    sol = M\rhs;
                end
                % Compute corrected values
                dalpha = sol(1:2*m);
                dt = rhop2-dalpha;
                dz = rhokkt1 - (z.*dalpha)./alpha;
                ds = rhokkt2 + (s.*dalpha)./t;
                dh = sol(end);
                
                % Compute lamdba factor
                comp = (e-1)*[alpha./dalpha t./dt];% s./ds z./dz];
                nu = min(1,min(comp(comp>=0)));
                if isempty(nu)
                    nu = 1;
                end
                
                
                % Compute new mu
                mu = (alpha'*z/(2*m))*((1.5-nu+e)/(10+nu))^2;
                
                % Update variables
                alpha = alpha + nu*dalpha;
                t = t + nu*dt;
                s = s + nu*ds;
                z = z + nu*dz;
                h = h + nu*dh;
                
                % Compute cancellation check
                gap = alpha'*z / abs(.5*alpha'*Q*alpha+c'*alpha+.5*alpha'*z);
                cnt = cnt+1;
                fprintf('Gap at iteration %d: %f, nu:%2.10f, mu:%f, used cholesky:%d\n',cnt,gap,nu,mu,cholflag);
            end
            if nu < eps
                disp('Stopped because nu got too small.');
            end
            % Return correct coefficients!
            ai = T*alpha;
            b = -h;
        end
    end
    
end

