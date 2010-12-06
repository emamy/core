classdef ScalarNuSVR < general.regression.BaseScalarSVR
    %SCALARSVR Scalar support vector regression.
    %
    %  Implementation details can be found in Daniel's Scratch
    %  Tex-Collection; it basically combines aspects from the books
    %  B. Schölkopf & A. Smola's "Learning with Kernels" and
    %  "Support Vector Machines" from I. Steinwart & A. Christman
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    properties
        % The weighting factor for epsilon against slack variables.
        % Value may range between 0 < nu < 1.
        %
        % Default: .4
        % See also: C
        nu = .4;
    end
    
    methods
        
        function this = ScalarNuSVR
            this = this@general.regression.BaseScalarSVR;
        end
        
        function [ai,b,svidx,epsi] = regress(this, fxi)
            %SCALAR_SVR Performs scalar support vector regression
            %
            % See also: KerMor
            
            % Compile quadratic program
            % Total number of samples
            m = size(this.K,1);
            % Ensure fxi is a column vector
            fxi = reshape(fxi,m,[]);
            
            % Check for zero function
            if all(abs(fxi) < eps)
                svidx = 1;
                ai = 0;
                b = 0;
                epsi = 0;
                warning('ScalarNuSvr:ZeroFun','All sample y_i values are less than eps. Assuming zero function.');
                return;
            end
            
            % Problem setup
            
            % Storage: alpha(1..m) = alpha_i, alpha(m+1..2m) = alpha_i*
            % T performs alpha_i* - alpha_i each
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            
            Q = T'*this.K*T;
            c = - T'*fxi;
            A = [ones(1,m) -ones(1,m); 
                 ones(1,2*m)];
            
            % Starting point
            %x0 = ones(2*m,1)*this.C/(2*m);
            x0 = [];
            
            % Bounds
            lbA = [0; 0];
            ubA = [0; this.C * this.nu];
            lb = zeros(2*m,1);
            ub = ones(2*m,1)*(this.C/m);
            
            % Call solver
            [p,d,info] = this.QPSolver.solve(Q,c,lb,ub,A,lbA,ubA,x0);
            
            % Convert results
            ai = T*p;
            b = -d(end-1);
            epsi = -d(end);
            svidx = find(abs(ai) >= this.AlphaMinValue);
            if isempty(svidx) && any(fxi ~= 0)
                error('No support vectors found. Problem unsolvable with current config?\nQuadprog exit flag:%d\nQuadprog out.message:%s',exitflag,out.message);
            end
            ai = ai(svidx);
            
            
%             %% Find "in-bound" SV's for eps/b computation
%             aiP = a(1:m); % \alpha_i^+
%             aiP=aiP(svidx);
%             aiN = a(m+1:end); % \alpha_i^-
%             aiN=aiN(svidx);
%             % Compile reduced versions of K and fxi
%             k = this.K(svidx,svidx);
%             svfxi = fxi(svidx);
%             
%             % Find the coefficients closest to C/2m
%             [v,Pidx] = min(abs(aiP-(this.C/(2*m))));
%             % Use the one closest to the upper bound if ai is not inside
%             % interval
%             useupper = true;
%             if aiP(Pidx) == 0 || aiP(Pidx) == this.C/m
%                 disp('Rare nuSVR case, no test written yet!');
%                 [v,Pidx] = min(svfxi' - ai'*k);
%                 useupper = false;
%             end
%             % Same for lower tube bounds
%             [v,Nidx] = min(abs(aiN-(this.C/(2*m))));
%             if aiN(Nidx) == 0 || aiN(Nidx) == this.C/m
%                 disp('Rare nuSVR case, no test written yet!');
%                 [v,Nidx] = min(ai'*k - svfxi');
%             end
%             
%             % Compute eps, b
%             epsi = (svfxi(Pidx)-svfxi(Nidx)-ai'*k(:,Pidx)+ai'*k(:,Nidx))/2;
%             %if epsi < 0
%                 %warning('some:id','nu-SVR failure! eps > 0 required.');
%                 disp(epsi);
%             %end
%             % Pick one of the indices to compute b; ideally one which had a
%             % coefficient inside the ]0,C/m[ interval
%             if useupper
%                 b = svfxi(Pidx) - ai'*k(:,Pidx) - epsi;
%             else
%                 b = svfxi(Nidx) - ai'*k(:,Nidx) + epsi;
%             end
            
            %poss = this.C/m-abs(ai) > this.AlphaMinValue;
            %bposs = [abs(ai(poss)'); svfxi(poss)' - ai'*k(:,poss) - sign(ai(poss))'*eps]
        end
        
        function set.nu(this, value)
            if ~isposrealscalar(value)
                error('nu must be a positive scalar');
            elseif value >= 1
                error('nu must be lower than one');
            end
            this.nu = value;
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            % Create new instance
            copy = general.regression.ScalarNuSVR;
            % Call superclass clone
            copy = clone@general.regression.BaseScalarSVR(this, copy);
            % Copy local props
            copy.nu = this.nu;
        end
    end
    
    methods(Static)
        function res = test_ScalarNuSVR
            % Performs a test of this class
            
            x = -5:.1:5;
            fx = sinc(x);
            %x = 1:10;
            %fx = ones(size(x))*5;
            
            svr = general.regression.ScalarNuSVR;
            svr.AlphaMinValue = 10*eps;
            svr.nu = .5;
            svr.C = 100000;
            %kernel = kernels.PolyKernel(2);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            
            figure(2);
            [ai,b,svidx,epsi] = svr.regress(fx);
            plot(x,fx,'r',x,[fx-epsi; fx+epsi],'r--');
            
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x) + b;
            
            fsvr = svfun(x);
            

            
            % Plot approximated function
            hold on;
            plot(x,fsvr,'b',x,[fsvr-epsi; fsvr+epsi],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');

            fdiff = abs(fsvr(svidx)-fx(svidx));
            res = isempty(find(fdiff ./ max(fx(svidx)) > 1e-7,1));
            if ~res
                errors = find(fdiff > epsi);
                plot(x(svidx(errors)),fx(svidx(errors)),'blackx','LineWidth',4);
            end
            
            tit = sprintf('nu = %f, b = %f, epsilon = %f, #SV=%d',svr.nu,b,epsi,length(svidx));
            title(tit);
            disp(tit);
            hold off;
            
            
        end
        
        function res = test_ScalarNuSVR_to_EpsSVR
            % Tests with random outliers
            
            x = -5:.1:5;
            
            r = rand(1,length(x));
            fx = sinc(x);
            fx(r<.2) = fx(r<.2)-.5;
            fx(r>.8) = fx(r>.8)+.5;
            
            %fx = sinc(x) + (rand(1,length(x))-.5)*.3;
            
            svr = general.regression.ScalarNuSVR;
            svr.nu = 0.3;
            svr.C = 10;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            [ai,b,svidx,epsi] = svr.regress(fx);
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x) + b;
            
            % Create eps-SVR and feed with computed epsilon
            esvr = general.regression.ScalarEpsSVR;
            esvr.eps = epsi;
            esvr.K = svr.K;
            esvr.C = svr.C;
            [eai,eb,esvidx,eeps] = esvr.regress(fx);
            esv = x(:,esvidx);
            esvfun = @(x)eai'*kernel.evaluate(esv,x) + eb;
            
            figure(1);
            xp = x(1,:);
            subplot(1,2,1);
            plot(xp,fx,'r');%,xp,[fx-eps; fx+eps],'r--');
            fsvr = svfun(x);
            hold on;
            
            % Plot approximated function
            plot(xp,fsvr,'b',xp,[fsvr-epsi; fsvr+epsi],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv(1,:),fx(svidx),'.r',xp(skipped),fx(skipped),'xr');
            title(sprintf('nu = %f, b = %f, epsilon = %f, #SV = %d, C = %d',svr.nu,b,epsi,length(svidx),svr.C));
            hold off;
            subplot(1,2,2);
            plot(xp,fx,'r');%,xp,[fx-eps; fx+eps],'r--');
            efsvr = esvfun(x);
            hold on;
            
            % Plot approximated function
            plot(xp,efsvr,'b',xp,[efsvr-eeps; efsvr+eeps],'b--');
            skipped = setdiff(1:length(x),esvidx);
            plot(esv(1,:),fx(esvidx),'.r',xp(skipped),fx(skipped),'xr');
            title(sprintf('b = %f, epsilon = %f, #SV = %d, C = %d',eb,eeps,length(esvidx),esvr.C));
            hold off;
            
            res = norm(fsvr-efsvr) < sqrt(eps);
        end
    end
end

