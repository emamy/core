classdef ScalarSVR < handle
    %SCALARSVR Scalar support vector regression.
    %
    %  See B. Schölkopf & A. Smola's "Learning with Kernels" for
    %  implementation details.
    %
    % @author Daniel Wirtz @date 11.03.2010
    
    properties
        % The kernel matrix to use.
        %
        % The reason why this is a property and not an argument is that
        % once a matrix is set multiple regressions for the same base
        % vector set can be performed easily.
        K;
        
        % The margin for the approximation-to-source distance.
        % Effectiveness also depends on C
        %
        % See also: C
        eps = 0.05;
        
        % The weighting of the slack variables. For larger C, the slack
        % variables are forced towards zero so that violations of the
        % eps-tube are getting penalized harder.
        % Theoretically, for C=Inf all fxi must be inside the eps-tube
        % around the original function.
        %
        % See also: eps
        C = 10;
        
        % Options for quadprog-solver
        QuadProgOpts = optimset('LargeScale','off','MaxIter',300,'Display','off');
        
        % Minimum value for any alpha to be considered a support vector
        % coefficient
        AlphaMinValue = 1e-10;
        
        % Maximum number of iterations for automatic b-computation in case
        % all alpha_i's are bounded by C/m
        bCompMaxIt = 100;
    end
    
    methods
        function [ai,b,svidx] = regress(this, fxi)
            %SCALAR_SVR Performs scalar support vector regression
            %
            % See also: globalconf
            
            % Compile quadratic program
            % Total number of samples
            m = size(this.K,1);
            % Ensure fxi is a column vector
            fxi = reshape(fxi,m,[]);
            
            % "Bugfix" for epsilons that are set larger than the range of
            % the function values to approximate.
            fxirange = max(fxi)-min(fxi);
            oldeps = [];
            if 2*this.eps >= fxirange
                oldeps = this.eps;
                this.eps = fxirange*.49;
                warning('ScalarSVR:epsTooLarge','eps is too large for data to approximate. Using eps=%1.5f temporarily.',this.eps);
            end
            
            % Storage: alpha(1..m) = alpha_i, alpha(m+1..2m) = alpha_i*
            % T performs alpha_i* - alpha_i each
            T = [-diag(ones(1,m)) diag(ones(1,m))];
            
            %% Problem setup
            prog.H = T'*this.K*T;
            prog.f = this.eps*ones(2*m,1) - T'*fxi;
            
            prog.Aeq = [ones(1,m) -ones(1,m)];
            prog.beq = 0;
            
            prog.lb = zeros(2*m,1);
            prog.ub = ones(2*m,1)*(this.C/m);
            
            prog.Aineq = [];
            prog.bineq = [];
            %prog.x0 = rand(2*m,1)*(this.C/m);
            
            prog.solver = 'quadprog';
            prog.options = this.QuadProgOpts;
            
            %% Solve quadratic problem
            %ai = quadprog(prog);
            % @TODO compile return messages for simulation summary/warnings
            [ai, fval, exitflag, out] = quadprog(prog);
            
            
            %% Extract support vectors
            % reduce ai from 2m to m vector
            % follow alpha_i* - alpha_i
            ai = T*ai;
            
            % Find support vectors
            svidx = find(abs(ai) >= this.AlphaMinValue);
            
            if isempty(svidx) && any(fxi ~= 0)
                error('No support vectors found. Problem unsolvable with current config?\nQuadprog exit flag:%d\nQuadprog out.message:%s',exitflag,out.message);
            end
            
            % check if b can be computed correctly
            if all(abs(ai(svidx)) - this.C/m < 1e-10)
                b = 0;
                
                % Save to-leave-out vectors
                skipped = setdiff(1:m,svidx);
                
                if ~isempty(skipped)
                    fxiskip = fxi(skipped)';
                    
                    ai = ai(svidx);
                    
                    oldb = Inf;
                    cnt = 0;
                    while abs(oldb-b) > 1e-5 && cnt < this.bCompMaxIt
                        % Compute approx fx values at xskipped
                        fsvrskip = ai'*this.K(svidx,skipped) + b;
                        if isempty(fsvrskip)
                            b = 4;
                        end
                        % Comput max and min differences
                        up  = max(fxiskip - fsvrskip + this.eps);
                        low = min(fxiskip - fsvrskip - this.eps);
                        
                        oldb = b;
                        b = b + (up + low) / 2;
                        cnt = cnt + 1;
                    end
                    if cnt == this.bCompMaxIt
                        warning('KerMor:svr:Ambiguous_offset',['All coefficients for SVR expansion are bounded.\n'...
                            'The offset b cannot be computed sufficiently precise after %d iterations, setting b=%f.'],...
                            this.bCompMaxIt,b);
                    end
                else
                    val = fxi' - ai'*this.K(svidx,:);
                    [tmp,i] = max(abs(val));
                    b = val(i);
                    fprintf('ScalarSVR warning: no skipped source vectors found.\nSetting b to max abs distance (=%f)\n',b);
                    %warning('KerMor:svr:noSkippedSupport',['No skipped source'...
                    %    'vectors found. Setting b to max abs distance (=%f).'],b);
                end
            else
                % Exact computation of b is possible due to the existance of an alpha_i
                % whose support vector lies on the border of the eps-tube.
                ai = ai(svidx);
                kmat = this.K(svidx,svidx);
                % Compute b
                [val,idx] = min(abs(abs(ai)-this.C/(2*m)));
                % get correct index regarding the complete input set
                index = svidx(idx(1));
                % compute b
                b = fxi(index) - ai' * kmat(:,idx(1)) + sign(m-index)*this.eps;
            end
            
            % Reset eps to the old value if changed at beginning.
            if ~isempty(oldeps)
                this.eps = oldeps;
            end
            
        end
        
        function set.K(this, value)
            % Make matrix symmetric (can be false due to rounding errors)
            this.K = .5*(value + value');
        end
    end
    
    methods(Static)
        function test_ScalarSVR
            % Performs a test of this class
            
            x = -5:.1:5;
            fx = sinc(x);
            %x = 1:10;
            %fx = ones(size(x))*5;
            
            svr = general.regression.ScalarSVR;
            svr.eps = .2;
            svr.C = 10;
            %kernel = kernels.PolyKernel(2);
            %kernel = kernels.LinearKernel;
            kernel = kernels.GaussKernel(1);
            svr.K = kernel.evaluate(x,x);
            
            figure(1);
            plot(x,fx,'r',x,[fx-svr.eps; fx+svr.eps],'r--');
            
            [ai,b,svidx] = svr.regress(fx);
            sv = x(:,svidx);
            svfun = @(x)ai'*kernel.evaluate(sv,x) + b;
            
            fsvr = svfun(x);
            
            hold on;
            
            % Plot approximated function
            plot(x,fsvr,'b',x,[fsvr-svr.eps; fsvr+svr.eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(sv,fx(svidx),'.r',x(skipped),fx(skipped),'xr');
            
            hold off;
        end
        
        function test2_CustomSVR
            % Performs a test of this class
            
            x = [0    0.3000    0.6000    0.9000;
                 0.8321    1.0756    1.3557    1.6539];
            fx = [0.7393    0.8799    0.9769    0.9966];
            
            svr = general.regression.ScalarSVR;
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
    end
    
end

