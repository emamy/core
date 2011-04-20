classdef ScalarEpsSVR < general.regression.BaseScalarSVR
    %SCALARSVR Scalar support vector regression.
    %
    %  Implementation details can be found in Daniel's Scratch
    %  Tex-Collection; it basically combines aspects from the books
    %  B. Schölkopf & A. Smola's "Learning with Kernels" and
    %  "Support Vector Machines" from I. Steinwart & A. Christman
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    properties
        
        % The margin for the approximation-to-source distance.
        % Effectiveness also depends on C
        %
        % See also: C
        eps = 0.05;
    end
    
    methods
        
        function this = ScalarEpsSVR
            % Default constructor. Calls superconstructor to initialize
            % properties
            this = this@general.regression.BaseScalarSVR;
        end
        
        function [ai,b,svidx,eps] = regress(this, fxi)
            % Performs regression using the IPOPT package
            
            m = size(this.K,1);
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            % Ensure fxi is a column vector
            fxi = reshape(fxi,m,[]);
            
            % Problem setup
            Q = T'*this.K*T;
            c = (this.eps*ones(2*m,1) - T'*fxi);
            A = [ones(1,m) -ones(1,m)];
            
            % Starting point
            %x0 = ones(2*m,1)*this.C/(2*m);
            x0 = [];
            
            % Bounds
            lbA = 0;
            ubA = 0;
            lb = zeros(2*m,1);
            ub = ones(2*m,1)*(this.C/m);
            
            % Call solver
            [p,d,info] = this.QPSolver.solve(Q,c,lb,ub,A,lbA,ubA,x0);
            
            % Convert results
            ai = T*p;
            b = -d(end);
            svidx = find(abs(ai) >= this.AlphaMinValue);
            ai = ai(svidx);
            eps = this.eps;
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            % Create new instance
            copy = general.regression.ScalarEpsSVR;
            % Call superclass clone
            copy = clone@general.regression.BaseScalarSVR(this, copy);
            % Copy local props
            copy.eps = this.eps;
        end
    end
    
    methods(Static)
        function res = test_ScalarEpsSVR
            % Performs a test of this class
            
            x = -5:.1:5;
            fx = sinc(x);
            %fx(30) = .5;
            %x = 1:10;
            %fx = ones(size(x))*5;
            
            svr = general.regression.ScalarEpsSVR;
            %svr.eps = 0.073648;
            svr.eps = .1;
            svr.C = 10;
            %svr.QPSolver.MaxIterations = 1000;
            %svr.QPSolver = solvers.qp.qpMatlab;
            %svr.QPSolver = solvers.qp.qpMosek;
            svr.QPSolver = solvers.qp.qpOASES;
            %kernel = kernels.PolyKernel(2);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(.7);
            svr.K = kernel.evaluate(x,x);
            
            figure;
            plot(x,fx,'r',x,[fx-svr.eps; fx+svr.eps],'r--');
            
            [ai,b,svidx,epsi] = svr.regress(fx);
            sv = x(:,svidx);
            svfun = @(x)ai'*(kernel.evaluate(x,sv)') + b;
            
            fsvr = svfun(x);
            
            fdiff = abs(fsvr(svidx)-fx(svidx));
            errors = find(fdiff  < .999*svr.eps);
            res = isempty(errors);
            
            % Plot approximated function
            hold on;
            plot(x,fsvr,'b',x,[fsvr-svr.eps; fsvr+svr.eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');
            
            if ~res
                plot(x(svidx(errors)),fx(svidx(errors)),'blackx','LineWidth',4);
            end
            
            tit = sprintf('#SV=%d, eps=%f, b=%f',length(svidx),epsi,b);
            title(tit);
            disp(tit);
            hold off;
        end
        
        function test2_CustomSVR
            % Performs a test of this class
            
            x = [0    0.3000    0.6000    0.9000;
                0.8321    1.0756    1.3557    1.6539];
            fx = [0.7393    0.8799    0.9769    0.9966];
            
            svr = general.regression.ScalarEpsSVR;
            svr.eps = 0.5;
            svr.C = 10;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            
            figure(1);
            xp = x(1,:);
            plot(xp,fx,'r',xp,[fx-svr.eps; fx+svr.eps],'r--');
            
            [ai,b,svidx] = svr.regress(fx);
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x) + b;
            
            fsvr = svfun(x);
            hold on;
            
            % Plot approximated function
            plot(xp,fsvr,'b',xp,[fsvr-svr.eps; fsvr+svr.eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv(1,:),fx(svidx),'.r',xp(skipped),fx(skipped),'xr');
            
            hold off;
        end
        
%         function res = test_SMOEpsSVR
%             
%             x = -5:.1:5;
%             fx = sinc(x);
%             ep = .1;
%             
%             % SVMToolbox-Part
%             kernel = rbf(ep);
%             C      = 1.0;
%             tutor  = smosvctutor;
%             net = train(svc, tutor, x', fx', C, kernel);
%             net = fixduplicates(net, x', fx');
%             net = strip(net);
%             sv = getsv(net);
%             ai = getw(net);
%             b = getbias(net);
%             fsvr = ai*evaluate(kernel,sv,x') + b;
%             
%             % Plotting
%             figure(1);
%             plot(x,fx,'r',x,[fx-ep; fx+ep],'r--');
%             
%             [dummy, skipped] = setdiff(int32(round(10*x)),int32(round(10*sv)));
%             svidx = setdiff(1:length(x),skipped);
%             
%             fdiff = abs(fsvr(svidx)-fx(svidx));
%             errors = find(fdiff  < .9999*ep);
%             res = isempty(errors);
%             
%             % Plot approximated function
%             hold on;
%             plot(x,fsvr,'b',x,[fsvr-ep; fsvr+ep],'b--');
%             
%             plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');
%             
%             if ~res
%                 plot(x(svidx(errors)),fx(svidx(errors)),'blackx','LineWidth',4);
%             end
%             
%             tit = sprintf('#SV=%d, eps=%f, b=%f',length(svidx),ep,b);
%             title(tit);
%             disp(tit);
%             hold off;
%         end
    end
    
end

