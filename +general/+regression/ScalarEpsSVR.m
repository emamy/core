classdef ScalarEpsSVR < general.regression.BaseQPSVR
%ScalarEpsSVR: Scalar support vector regression.
%
%  Implementation details can be found in Daniel's Scratch
%  Tex-Collection; it basically combines aspects from the books
%  B. Sch�lkopf & A. Smola's "Learning with Kernels" and
%  "Support Vector Machines" from I. Steinwart & A. Christman
%
% See also: ScalarNuSVR
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
% @change{0,4,sa,2011-05-06} Implemented Setter for the property eps
%
% @change{0,4,dw,2011-05-03} Removed the `b` offset terms from the SVM formulation (i.e.
% removing the equality constraint for the coefficient vectors) as kernel expansions no longer
% use offset terms.
    
    properties(SetObservable)
        % The margin for the approximation-to-source distance.
        % Effectiveness also depends on C
        %
        % @propclass{critical} Determines the precision of the
        % approximation.
        %
        % @default 0.05
        %
        % See also: C
        Eps = 0.05;
    end
    
    methods
        function this = ScalarEpsSVR
            % Default constructor. Calls superconstructor to initialize
            % properties
            this = this@general.regression.BaseQPSVR;
            this.registerProps('eps');
        end
        
        function ai = regress(this, fxi, ainit)
            % Performs epsilon-SVR regression.
            %
            % Parameters:
            % fxi: The `f(x_i)` values to regress given the kernel matrix `K` computed from
            % `\Phi(x_i,x_j)`.
            % ainit: [Optional] An initial value set for the coefficients.
            %
            % Return values:
            % ai: The kernel expansion coefficients `\alpha_i`.
            m = size(this.K,1);
            T = [diag(ones(1,m)) -diag(ones(1,m))];
            
            % Problem setup
            Q = T'*this.K*T;
            c = (this.Eps*ones(2*m,1) - (fxi*T)');
            
            % Starting point
            if nargin < 3
                %ainit = ones(2*m,1)*this.C/(2*m);
                ainit = [];
            end
            
            % Bounds
            lb = zeros(2*m,1);
            ub = ones(2*m,1)*(this.C/m);
            
            % Call solver
            p = this.solve(Q,c,lb,ub,[],[],[],ainit);
            
            % Convert results
            ai = T*p;
        end
        
        function set.Eps(this, value)
            if ~isposrealscalar(value)
                error('Value must be a positive real scalar');
            end
            this.Eps = value;
        end
    end
    
    methods(Sealed)
        function copy = clone(this)
            % Create new instance
            copy = general.regression.ScalarEpsSVR;
            % Call superclass clone
            copy = clone@general.regression.BaseQPSVR(this, copy);
            % Copy local props
            copy.Eps = this.Eps;
        end
    end
    
    methods(Static)
        function res = test_ScalarEpsSVR
            % Performs a test of this class
            
            x = -5:.1:5;
            fx = sinc(x)-1;
            %fx(30) = .5;
            %x = 1:10;
            %fx = ones(size(x))*5;
            
            svr = general.regression.ScalarEpsSVR;
            %svr.eps = 0.073648;
            svr.Eps = .1;
            svr.Lambda = 1/20;
            %svr.MaxIterations = 1000;
            
            %kernel = kernels.PolyKernel(2);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(.5);
            svr.K = kernel.evaluate(x,x);
            
            figure;
            plot(x,fx,'r',x,[fx-svr.Eps; fx+svr.Eps],'r--');
            
            [ai, svidx] = svr.computeKernelCoefficients(fx, []);
            sv = x(:,svidx);
            svfun = @(x)ai'*(kernel.evaluate(x,sv)');
            
            fsvr = svfun(x);
            
            fdiff = abs(fsvr(svidx)-fx(svidx));
            errors = find(fdiff  < .999*svr.Eps);
            res = isempty(errors);
            
            % Plot approximated function
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
        
        function test2_CustomSVR
            % Performs a test of this class
            
            x = [0    0.3000    0.6000    0.9000;
                0.8321    1.0756    1.3557    1.6539];
            fx = [0.7393    0.8799    0.9769    0.9966];
            
            svr = general.regression.ScalarEpsSVR;
            svr.Eps = 0.5;
            svr.C = 10;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            
            figure(1);
            xp = x(1,:);
            plot(xp,fx,'r',xp,[fx-svr.Eps; fx+svr.Eps],'r--');
            
            [ai,svidx] = svr.regress(fx);
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x);
            
            fsvr = svfun(x);
            hold on;
            
            % Plot approximated function
            plot(xp,fsvr,'b',xp,[fsvr-svr.Eps; fsvr+svr.Eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv(1,:),fx(svidx),'.r',xp(skipped),fx(skipped),'xr');
            
            hold off;
        end
    end
end

