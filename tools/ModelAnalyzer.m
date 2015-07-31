classdef ModelAnalyzer < handle;
%ModelAnalyzer: Analysis tools for reduced models and approximations
%
% @author Daniel Wirtz @date 2011-11-17
%
% @change{0,6,dw,2011-11-17} Moved this class to the +tools package from visual.
%    
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        
        function [params, errs, details] = getRedErrForRandomParamSamples(this, num, seed, in)
            % Computes the simulation errors (output) for 'num' random model parameters.
            %
            % Parameters:
            % num: The number of random parameters to get the error for. @type integer @default
            % 50
            % seed: The seed for the random number generator. @type integer @default
            % 'round(cputime*100)'
            % in: The number of the input to use. Leave unset for no model inputs.
            %
            % Return values:
            % params: The parameters generated. @type matrix<double>
            % errs: A `4\times num` matrix containing the Linf-L2 absolute and relative error
            % in rows 1,2 and Linf-Linf absolute and relative errors in rows 3-4.
            %
            % See also:
            % ModelAnalyzer.getRedErrForParamSamples
            if nargin < 4
                in = [];
                if nargin < 3
                    seed = round(cputime*100);
                    if nargin < 2
                        num = 50;
                    end
                end
            end
            params = this.rm.FullModel.getRandomParam(num, seed);
            [errs, details] = getRedErrForParamSamples(this, params, in);
        end
        
        function [errs, details] = getRedErrForParamSamples(this, params, in)
            % Computes the simulation errors (output) for the given model parameters.
            %
            % Parameters:
            % params: The parameters to sample and compute
            % the error for. If not given, the param samples for reduced model computation are
            % used.  @type integer @default FullModel.Data.ParamSamples
            % in: The input index to use for the simulations. @type integer
            % @default []
            %
            % Return values:
            % errs: A `5\times n` matrix containing the Linf-L2 absolute and relative error in
            % rows 1,2 and Linf-Linf absolute and relative errors in rows 3-4. Row 5 contains
            % the error estimator efficiency if enabled.
            fm = this.rm.FullModel;
            if nargin < 3
                in = [];
                if nargin < 2
                    params = fm.Data.ParamSamples;
                end
            end
            num = size(params,2);
            details = struct;
            details.params = params;
            details.L2Errs = zeros(num,length(fm.Times));
            details.L2Fullsol = details.L2Errs;
            details.L2OutErrs = details.L2Errs;
            details.L2OutFullsol = details.L2Errs;
            details.Linferrs = details.L2Errs;
            details.Linftruesolnorm = details.L2Errs;
            details.Times = zeros(num,2); % computation times full/reduced
            haveerrest = ~isempty(this.rm.ErrorEstimator) && this.rm.ErrorEstimator.Enabled;
            if haveerrest
                details.Estimates = details.L2Errs;
                details.OutEstimates = details.L2Errs;
                errs = zeros(8,num);
            else
                errs = zeros(6,num);
            end
            pi = ProcessIndicator('Computing reduction error details for %d param samples',...
                num,false,num);
            for pidx = 1:num
                mu = params(:,pidx);
                [~, y, tf, x] = fm.simulate(mu, in);
                [~, yr, tr, xr] = this.rm.simulate(mu, in);
                details.Times(pidx,:) = [tf tr];
                details.L2Errs(pidx,:) = Norm.L2(x-this.rm.V*xr);
                details.L2Fullsol(pidx,:) = Norm.L2(x);
                details.L2OutErrs(pidx,:) = Norm.L2(y-yr);
                details.L2OutFullsol(pidx,:) = Norm.L2(y);
                details.Linferrs(pidx,:) = Norm.Linf(x-this.rm.V*xr);
                details.Linftruesolnorm(pidx,:) = Norm.Linf(x);
                errs(1,pidx) = max(details.L2Errs(pidx,:)); %linf l2 err
                errs(2,pidx) = errs(1,pidx) ./ max(details.L2Fullsol(pidx,:)); % rel
                errs(3,pidx) = max(details.L2OutErrs(pidx,:)); %linf l2 err - output
                errs(4,pidx) = errs(3,pidx) ./ max(details.L2OutFullsol(pidx,:)); % rel - output
                errs(5,pidx) = max(details.Linferrs(pidx,:)); %linf linf err
                errs(6,pidx) = errs(5,pidx) ./ max(details.Linftruesolnorm(pidx,:)); % rel
                if haveerrest
                    % Effectivity (w.r.t. L2 error)
                    % State space
                    details.Estimates(pidx,:) = this.rm.ErrorEstimator.StateError;
                    errs(5,pidx) = details.Estimates(pidx,end)/details.L2Errs(pidx,end);
                    % Output
                    details.OutEstimates(pidx,:) = this.rm.ErrorEstimator.OutputError;
                    errs(6,pidx) = details.OutEstimates(pidx,end)/details.L2OutErrs(pidx,end);
                end
                pi.step;
            end
            pi.stop;
        end
        
        function plotRedErrForParamSamples(this, errs, pm)
           % @\todo implement, so far only plotting in DEIM testing with multiple DEIM m Orders.
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
                str = sprintf('%s, mu=[%s]',str,Utils.implode(mu,', ','%2.3f'));
            end
            if ~isempty(inputidx)
                str = sprintf('%s, u_%d',str,inputidx);
            end
            t = PrintTable('Computation times %s:',str);
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
            %fprintf('||y(t_i)||_2: %s',Utils.implode(l2,', ','%2.3f'));
            t = PrintTable('Error comparison for %s:',str);
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
            if nargout < 1
                t.display;
            end
            
            %% Plotting
            pm = PlotManager(false, 2, 2);
            fm.plot(ti, y, pm.nextPlot('full_sim',sprintf('Full simulation\n%s',str)));
            fm.plot(ti, yr, pm.nextPlot('red_sim',sprintf('Reduced simulation\n%s',str)));
            fm.plot(ti, log10(abs(y-yr)), pm.nextPlot('abs_err',sprintf('Absolute error (log scale)\n%s',str)));
            hlp = abs(y);
            fm.plot(ti, log10(abs(y-yr)./hlp), pm.nextPlot('rel_err',sprintf('Relative error (log scale)\n%s',str)));
            pm.done;
            if nargout < 2
                pm.LeaveOpen = true;
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
                if ~isempty(this.rm.V)
                    x = this.rm.V*(this.rm.W'*x);
                end
                mu = repmat(mu,1,numel(t));
                fx = fm.Approx.evaluateMulti(x,t,mu);
                afx = fm.System.f.evaluateMulti(x,t,mu);
                el2 = Norm.L2(fx-afx);
                elinf = Norm.Linf(fx-afx);
            else
                error('The approximation error can only be computed for models with an approx.BaseApprox instance present.');
            end
        end
        
        function [el2, elinf, pm] = getATDError(this, pm)
            % Computes the approximation training error on the full
            % system's trajectory for given mu and inputidx.
            
            rm = this.rm;
            fm = rm.FullModel;
            
            if ~isempty(fm.Approx)
                if nargin < 2
                    pm = PlotManager(false,2,2);
                end
                % Full approx
                atd = fm.Data.ApproxTrainData;
                afx =  fm.Approx.evaluateMulti(atd.xi.toMemoryMatrix,atd.ti,atd.mui);
                s = 1:size(atd.xi,2);
                fx = atd.fxi.toMemoryMatrix;
                nofx = Norm.L2(fx);
                el2 = Norm.L2(fx-afx);
                elinf = Norm.Linf(fx-afx);
                h = pm.nextPlot('abs_l2','Absolute L2 error','x_i','L2');
                LogPlot.cleverPlot(h,s,el2);
                h = pm.nextPlot('rel_l2','Relative L2 error','x_i','L2');
                LogPlot.cleverPlot(h,s,el2./nofx);
                
                if ~isempty(rm.V)
                    rd = size(rm.V);
                    % Projected variant
                    apfx =  this.rm.System.f.evaluateMulti(rm.W'*atd.xi,atd.ti,atd.mui);
                    nofx = Norm.L2(afx);
                    pel2 = Norm.L2(afx-rm.V*apfx);
                    h = pm.nextPlot('abs_proj_l2',...
                        sprintf('Absolute L2 error projected vs full approx, %d/%d dims',rd(2),rd(1)),...
                        'x_i','L2');
                    LogPlot.cleverPlot(h,s,pel2);
                    h = pm.nextPlot('rel_proj_l2',...
                        sprintf('Relative L2 error projected vs full approx, %d/%d dims',rd(2),rd(1)),...
                        'x_i','L2');
                    LogPlot.cleverPlot(h,s,pel2./nofx);
                end
                
                pm.done;
                if nargout < 3
                    pm.LeaveOpen = true;
                end
            else
                error('The approximation error can only be computed for models with an approx.BaseApprox instance present.');
            end
            
        end
        
        function pm = analyzeError(this, mu, inputidx, pm)
            if isempty(this.rm.ErrorEstimator)
                error('Error analysis only available for models with error estimator');
            end
            rmodel = this.rm;
            if nargin < 4
                pm = PlotManager(false, 2, 3);
                if nargin < 3
                    inputidx = rmodel.DefaultInput;
                    if nargin < 2
                        mu = rmodel.DefaultMu;
                    end
                end
            end
            
            
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
            plot(h,t,Utils.preparePlainPlot(x));
            axis([t(1) t(end) ymin ymax]);
            
            if this.UseOutput
                h = pm.nextPlot('outputs',...
                    sprintf('Full+reduced system outputs with error bounds\n(r=%d,self time:%.3f, time with err est:%.3f)',...
                    size(rmodel.V,2),timer_noerr,timer),...
                    'time','y(t) / y^r(t)');
                plot(h, t,Utils.preparePlainPlot(x),'b',...
                    t, Utils.preparePlainPlot(xr),'r',...
                    t, Utils.preparePlainPlot(xrmin),...
                    'r--',t,Utils.preparePlainPlot(xrplus),'r--');
                legend('Full system','Reduced system','Lower bound','Upper bound');
            else
                h = pm.nextPlot('redsys',sprintf('Reduced system (r=%d,self time:%.3f,\ntime with err est:%.3f)',size(rmodel.V,2),timer_noerr,timer),...
                    'time','x^r(t)');
                plot(h, t, Utils.preparePlainPlot(xr),'r');
            end
            axis([t(1) t(end) ymin ymax]);
            
            %% Absolute value plot
            h = pm.nextPlot('norms','The state variable norms','time','norms');
            plot(h, t, xnorm,'b',t,xrnorm,'r',t,xrnorm-est,'r--',t,xrnorm+est,'r--');
            le = legend('Full system','Reduced system','Lower bound','Upper bound');
            set(le,'Location','NorthWest');
            
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
        
        function pm = plotReductionOverview(this, pm)
            % @todo create interface with plot method -> current config gives automatic
            % reduction overview in each component
            if nargin < 2
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            [~,m,~,l] = FindInstance(this.rm.FullModel,'IReductionSummaryPlotProvider');
            for k=1:length(m)
                objs = m{k};
                nobj = length(objs);
                extra = '';
                for i=1:nobj
                    if nobj > 0
                        extra = sprintf('-%d',i);
                    end
                    context = sprintf('%s%s (%s)',l{k},extra,class(objs(i)));
                    objs(i).plotSummary(pm,context);
                end
            end
            pm.done;
        end
    end
end

