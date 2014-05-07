classdef ScalarNuSVR < general.regression.BaseQPSVR
%SCALARSVR Scalar support vector regression.
%
%  Implementation details can be found in Daniel's Scratch
%  Tex-Collection; it basically combines aspects from the books
%  B. Sch�lkopf & A. Smola's "Learning with Kernels" (p.260ff) and
%  "Support Vector Machines" from I. Steinwart & A. Christman
%
% See also: ScalarEpsSVR
%
% @author Daniel Wirtz @date 11.03.2010
%
% @change{0,7,dw,2013-01-23} Removed QP solvers and introduced intermediate class BaseQPSVR to
% include qp stuff.
%
% @change{0,5,dw,2011-09-09} Moved the QPSolver property to this class.
%
% @change{0,4,dw,2011-06-01} Fitted the interface to the changed one in
% general.regression.BaseScalarSVR
%
% @change{0,4,dw,2011-05-04} Removed offset `b` terms from computations.
    
    properties(SetObservable)
        % The weighting factor for epsilon against slack variables.
        % Value may range between 0 < nu < 1.
        %
        % The `\nu` SVR parameter determines the resulting `\epsilon` value. `\nu` is an upper bound
        % on the fraction of errors and a lower bound on the fraction of support vectors.
        %
        % @propclass{critical} Too low `\nu` parameters might result in too few support vectors and thus
        % large `epsilon`. Is closely connected to the C property.
        %
        % @default .4
        %
        % See also: C
        nu = .4;
    end
    
    properties(Transient, SetAccess=private)
        LastEpsilon;
    end
    
    methods
        
        function this = ScalarNuSVR
            this = this@general.regression.BaseQPSVR;
            this.registerProps('nu');
        end
        
        function [ai, sf] = regress(this, fxi, ainit)
            % Performs scalar nu-support vector regression
            %
            % Parameters:
            % fxi: The `f(x_i)` values to regress given the kernel matrix `K` computed from
            % `\K(x_i,x_j)`.
            % ainit: [Optional] An initial value set for the coefficients.
            %
            % Return values:
            % ai: The kernel expansion coefficients `\alpha_i`.
            
            %% Compile quadratic program
            
            % Total number of samples
            m = size(this.K,1);
            
            % Storage: alpha(1..m) = alpha_i, alpha(m+1..2m) = alpha_i*
            % T performs alpha_i* - alpha_i each
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            
            Q = T'*this.K*T;
            Q = (Q+Q')/2;
            c = -(fxi*T)';
            A = ones(1,2*m);
            
            % Starting point
            if nargin < 3 || isempty(ainit)
                ainit = ones(m,1)*this.C/m;
                %ainit = [];
            end
            
            % Bounds
            lbA = 0;
            ubA = this.C * this.nu;
            lb = zeros(2*m,1);
            ub = ones(2*m,1)*(this.C/m);
            
            %% Call solver
            [p,d,info] = this.solve(Q,c,lb,ub,A,lbA,ubA,T'*ainit);
            
            %% Convert results
            ai = T*p;
            this.LastEpsilon = d(end);
            
            %| @todo return correct stop flag
            sf = StopFlag.SUCCESS;
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
            svr.nu = .1;
            svr.Lambda = 1/100;
            %kernel = kernels.PolyKernel(2);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            
            figure(2);
            [ai, svidx] = svr.computeKernelCoefficients(fx, []);
            epsi = svr.LastEpsilon;
            plot(x,fx,'r',x,[fx-epsi; fx+epsi],'r--');
            
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x);
            
            fsvr = svfun(x);
            
            % Plot approximated function
            hold on;
            plot(x,fsvr,'b',x,[fsvr-epsi; fsvr+epsi],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');

            fdiff = abs(fsvr(svidx)-fx(svidx));
            errors = find(fdiff  < .999*epsi);
            res = isempty(errors);
            if ~res
                plot(x(svidx(errors)),fx(svidx(errors)),'blackx','LineWidth',4);
            end
            
            tit = sprintf('nu = %f, epsilon = %f, #SV=%d',svr.nu,epsi,length(svidx));
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
            svr.Lambda = 1/20;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            [ai,svidx] = svr.computeKernelCoefficients(fx,[]);
            epsi = svr.LastEpsilon;
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x);
            
            % Create eps-SVR and feed with computed epsilon
            esvr = general.regression.ScalarEpsSVR;
            esvr.Eps = epsi;
            esvr.K = svr.K;
            esvr.Lambda = svr.Lambda;
            [eai,esvidx] = esvr.computeKernelCoefficients(fx,[]);
            esv = x(:,esvidx);
            esvfun = @(x)eai'*kernel.evaluate(esv,x);
            
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
            title(sprintf('nu = %f, epsilon = %f, #SV = %d, C = %d',svr.nu,epsi,length(svidx),svr.C));
            hold off;
            subplot(1,2,2);
            plot(xp,fx,'r');%,xp,[fx-eps; fx+eps],'r--');
            efsvr = esvfun(x);
            hold on;
            
            % Plot approximated function
            plot(xp,efsvr,'b',xp,[efsvr-epsi; efsvr+epsi],'b--');
            skipped = setdiff(1:length(x),esvidx);
            plot(esv(1,:),fx(esvidx),'.r',xp(skipped),fx(skipped),'xr');
            title(sprintf('epsilon = %f, #SV = %d, C = %d',epsi,length(esvidx),esvr.C));
            hold off;
            
            res = norm(fsvr-efsvr) < sqrt(eps);
        end
    end
end

