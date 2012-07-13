classdef ModelAnalyzer < handle;
%ModelAnalyzer: Analysis tools for reduced models and approximations
%
% @author Daniel Wirtz @date 2011-11-17
%
% @change{0,6,dw,2011-11-17} Moved this class to the +tools package from visual.
%    
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing 
%
% @todo think of expressive names for methods
    
    properties
        SingleFigures = false;
        
        UseOutput = true;
    end
    
    properties(Access=private)
        rm;
    end
    
    methods
        
        function this = ModelAnalyzer(rmodel)
           if ~isa(rmodel,'models.ReducedModel')
               error('rmodel must be a models.ReducedModel subclass.');
           end
           this.rm = rmodel; 
        end
        
        function errs = getRedErrForParamSamples(this, in)
            % Computes the reduction error for all parameter samples in the full model's
            % ModelData.
            if nargin < 2
                in = [];
            end
            fm = this.rm.FullModel;
            errs = zeros(2,fm.Data.SampleCount);
            for pidx = 1:fm.Data.SampleCount
                mu = fm.Data.ParamSamples(:,pidx);
                y = fm.Data.getTrajectory(mu,in);
                [~, yr] = this.rm.simulate(mu,in);
                errs(1,pidx) = max(Norm.L2(yr-y)); %linf l2 err
                errs(2,pidx) = max(Norm.Linf(yr-y)); %linf linf err
            end
        end
        
        function errs = getRedErrForRandomParamSamples(this, num, in)
            % Computes the simulation errors (output) for 'num' random model parameters.
            %
            % Parameters:
            % num: The number `n` of random parameters to sample and compute
            % the error for. @type integer
            % in: The input index to use for the simulations. @type integer
            % @default []
            %
            % Return values:
            % errs: A `4\times n` matrix containing the Linf-L2 absolute
            % and relative error in rows 1,2 and Linf-Linf absolute and
            % relative errors in rows 3-4.
            if nargin < 3
                in = [];
            end
            fm = this.rm.FullModel;
            errs = zeros(4,num);
            for pidx = 1:num
                mu = fm.getRandomParam;
                [~, y] = fm.simulate(mu, in);
                [~, yr] = this.rm.simulate(mu, in);
                errs(1,pidx) = max(sqrt(sum((yr-y).^2))); %linf l2 err
                errs(2,pidx) = errs(1,pidx) / max(sqrt(sum(y.^2))); % rel
                errs(3,pidx) = max(max(abs(yr-y),[],1)); %linf linf err
                errs(4,pidx) = errs(3,pidx) / max(max(abs(y),[],1)); % rel
            end
        end
        
        function [t, pm] = compareRedFull(this, mu, inputidx)
            % Compares the solutions of the reduced model and the associated full model by
            % calling the BaseModel.plot method for both solutions and again for the
            % difference. Also some information of `l^2` and `l^\infty` errors are printed.
            %
            % Parameters:
            % mu: The concrete mu parameter sample to simulate for.
            % inputidx: The index of the input function to use.
            %
            % Return values:
            % t: The PrintTable instance
            if nargin < 3
                inputidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            fm = this.rm.FullModel;
            [~, y, ftime] = fm.simulate(mu,inputidx);
            [ti,yr, rtime] = this.rm.simulate(mu,inputidx);
            %% Text output
            str = sprintf('%s',fm.Name);
            if ~isempty(mu)
                str = sprintf('%s, mu=[%s]',str,general.Utils.implode(mu,', ','%2.3f'));
            end
            if ~isempty(inputidx)
                str = sprintf('%s, u_%d',str,inputidx);
            end
            t = PrintTable('Computation times %s:\n',str);
            t.HasRowHeader = true;
            t.addRow('Full model',ftime,{'%2.4fs'});
            t.addRow('Reduced model',rtime,{'%2.4fs'});
            t.addRow('Speedup',ftime/rtime,{'x%2.4f'});
            t.display;
            
            % L^2 errors
            l2 = Norm.L2(yr-y);
            lil2 = max(l2);
            l2l2 = Norm.L2(l2');
            meanl2 = mean(l2);
            l2relyl2 = l2 ./ Norm.L2(y);
            l2l2relyl2 = Norm.L2(l2relyl2');
            lil2relyl2 = max(l2relyl2);
            meanrell2 = mean(l2relyl2);
            %fprintf('||y(t_i)||_2: %s',general.Utils.implode(l2,', ','%2.3f'));
            t = PrintTable('Error comparison for %s:\n',str);
            t.addRow('L2 time and space error','L^2(||y(t) - yr(t)||_2,[0,T])',l2l2);
            t.addRow('Linf time and L2 space error','L^inf(||y(t) - yr(t)||_2,[0,T])',lil2);
            t.addRow('Relative L2 time and space error','L^2(||(y(t) - yr(t)) / y(t)||_2,[0,T])',l2l2relyl2);
            t.addRow('Relative Linf time and L2 space error', 'L^inf(||(y(t) - yr(t)) / y(t)||_2,[0,T])',lil2relyl2);
            t.addRow('Mean L2 error','Mean(||y(t) - yr(t)||_2,[0,T])',meanl2);
            t.addRow('Mean relative L2 error','Mean(||(y(t) - yr(t)) / y(t)||_2,[0,T])',meanrell2);
            
            % L^inf errors
            li = Norm.Linf(yr-y);
            lili = max(li);
            l2li = Norm.L2(li');
            meanli = mean(li);
            lirelyli = li ./ Norm.Linf(y);
            l2lirelyli = Norm.L2(lirelyli');
            lilirelyli = max(lirelyli);
            meanrelli = mean(lirelyli);
            t.addRow('Linf time and space error','L^inf(||y(t) - yr(t)||_inf,[0,T])',lili);
            t.addRow('L2 time and Linf space error','L^2(||y(t) - yr(t)||_inf,[0,T])',l2li);
            t.addRow('Relative Linf time and space error','L^inf(||(y(t) - yr(t)) / y(t)||_inf,[0,T])',lilirelyli);
            t.addRow('Relative L2 time and Linf space error','L^2(||(y(t) - yr(t)) / y(t)||_inf,[0,T])',l2lirelyli);
            t.addRow('Mean Linf error','Mean(||y(t) - yr(t)||_inf,[0,T])',meanli);
            t.addRow('Mean relative Linf error','Mean(||(y(t) - yr(t)) / y(t)||_inf,[0,T])',meanrelli);
            t.display;
            
            %% Plotting
            if nargout > 1
                pm = tools.PlotManager(false, 2, 2);
                fm.plot(ti, y, pm.nextPlot('full_sim',sprintf('Full simulation\n%s',str)));
                fm.plot(ti, yr, pm.nextPlot('red_sim',sprintf('Reduced simulation\n%s',str)));
                fm.plot(ti, log10(abs(y-yr)), pm.nextPlot('abs_err',sprintf('Absolute error (log scale)\n%s',str)));
                hlp = abs(y);
%                 if any(hlp(:) == 0)
%                     hlp2 = hlp;
%                     hlp2(hlp==0) = [];
%                     ep = min(hlp2(:))^2;
%                     hlp(hlp==0) = ep;
%                 end
                fm.plot(ti, log10(abs(y-yr)./hlp), pm.nextPlot('rel_err',sprintf('Relative error (log scale)\n%s',str)));
                pm.done;
            end
        end
        
        function [el2, elinf, t, x, fx, afx] = getTrajApproxError(this, mu, inputidx)
            % Computes the approximation training error on the full
            % system's trajectory for given mu and inputidx.
            fm = this.rm.FullModel;
            if ~isempty(fm.Approx)
                if nargin == 2
                    inputidx = [];
                end
                % Computes trajectory (including cache-lookup in
                % BaseFullModel!)
                [t,x] = fm.computeTrajectory(mu, inputidx);
                if ~isempty(fm.Data.V)
                    x = fm.Data.V*(fm.Data.W'*x);
                end
                mu = repmat(mu,1,numel(t));
                fx = fm.Approx.evaluate(x,t,mu);
                afx = fm.System.f.evaluate(x,t,mu);
                el2 = Norm.L2(fx-afx);
                elinf = Norm.Linf(fx-afx);
            else
                error('The approximation error can only be computed for models with an approx.BaseApprox instance present.');
            end
        end
        
        function pm = analyzeError(this, mu, inputidx, pm)
            if nargin < 4
                pm = tools.PlotManager(false, 2, 3);
                if nargin < 3
                    inputidx = [];
                    if nargin < 2
                        mu = [];
                    end
                end
            end
            rmodel = this.rm;
            
            %% Initial computations
            [~, y, time, x] = rmodel.FullModel.simulate(mu, inputidx);
            rmodel.ErrorEstimator.Enabled = false;
            tic; rmodel.simulate(mu, inputidx); timer_noerr = toc;
            rmodel.ErrorEstimator.Enabled = true;
            [t, yr, timer, xr] = rmodel.simulate(mu, inputidx);
            
            if this.UseOutput
                x = y;
                xr = yr;
                est = rmodel.ErrorEstimator.OutputError;
            else
                if ~isempty(rmodel.V)
                    xr = rmodel.V*xr;
                end
                est = rmodel.ErrorEstimator.StateError;
            end
            e = Norm.L2(x - xr);
            xnorm = Norm.L2(x);
            erel = e./xnorm;
            estrel = est./xnorm;
            xrnorm = Norm.L2(xr);
            erelr = e./xrnorm;
            estrelr = est./xrnorm;
            
            %% System plot
            xrmin = xr-repmat(est,size(xr,1),1); xrplus = xr+repmat(est,size(xr,1),1);
            ymax = max([max(x(:)) max(xr(:)) max(xrmin(:)) max(xrplus(:))]);
            ymin = min([min(x(:)) min(xr(:)) min(xrmin(:)) min(xrplus(:))]);
            
            if this.UseOutput
                ti = sprintf('The full system''s output (m=%d,time=%.3f)',size(x,1),time);
            else
                ti = sprintf('The full system (d=%d,time=%.3f)',size(x,1),time);
            end
            h = pm.nextPlot('fullsys',ti,'time','error');
            plot(h,t,general.Utils.preparePlainPlot(x));
            axis([t(1) t(end) ymin ymax]);
            
            if this.UseOutput
                h = pm.nextPlot('outputs',...
                    sprintf('Full+reduced system outputs with error bounds\n(r=%d,self time:%.3f, time with err est:%.3f)',...
                    size(rmodel.V,2),timer_noerr,timer),...
                    'time','y(t) / y^r(t)');
                plot(h, t,general.Utils.preparePlainPlot(x),'b',...
                    t, general.Utils.preparePlainPlot(xr),'r',...
                    t, general.Utils.preparePlainPlot(xrmin),...
                    'r--',t,general.Utils.preparePlainPlot(xrplus),'r--');
                legend('Full system','Reduced system','Lower bound','Upper bound');
            else
                h = pm.nextPlot('redsys',sprintf('Reduced system (r=%d,self time:%.3f,\ntime with err est:%.3f)',size(rmodel.V,2),timer_noerr,timer),...
                    'time','x^r(t)');
                plot(h, t, general.Utils.preparePlainPlot(xr),'r');
            end
            axis([t(1) t(end) ymin ymax]);
            
            %% Absolute value plot
            h = pm.nextPlot('norms','The state variable norms','time','norms');
            plot(h, t, xnorm,'b',t,xrnorm,'r',t,xrnorm-est,'r--',t,xrnorm+est,'r--');
            legend('Full system','Reduced system','Lower bound','Upper bound');
            
            % Error plots
            h = pm.nextPlot('abserrors',sprintf('The state variable absolute errors.\nmean(e)=%g, mean(est)=%g',mean(e),mean(est)),...
                'time','errors');
            semilogy(h,t,est,'b',t,e,'r');%,t,abs(e-est),'g');
            legend('Estimated error','True error');%,'Location','Best');
            
            % Relative Error plots
            h = pm.nextPlot('relerrors',sprintf(['The state variable relative errors (comp. to '...
                'full solution)\nmean(e_{rel})=%g, mean(est_{rel})=%g'],...
                mean(erel),mean(estrel)),...
                'time','errors');
            semilogy(h,t,estrel,'b',t,erel,'r');
            legend('Estimated error','True error');%,'Location','Best');
            
            h = pm.nextPlot('relerrors_red',sprintf(['The state variable relative errors (comp. to '...
                'reduced solution)\nmean(ered_{rel})=%g, mean(estred_{rel})=%g'],...
                mean(erelr),mean(estrelr)),...
                'time','errors');
            semilogy(h,t,erelr,'r',t,estrelr,'b');
            legend('True error','Estimated error');%,'Location','Best');
            if nargin < 4
                pm.done;
            end
        end
    end
end

