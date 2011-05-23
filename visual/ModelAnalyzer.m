classdef ModelAnalyzer < handle;
    %MODELANALYZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SingleFigures = false;
        
        UseOutput = true;
    end
    
    methods
        
        function this = ModelAnalyzer(rmodel, varargin)
            if nargin > 0
                this.analyze(rmodel,varargin{:});
            end
        end
        
        function analyze(this, rmodel, mu, inputidx)
            if nargin < 4
                inputidx = 1;
                if nargin < 3
                    mu = [];
                end
            end
            
            %% Initial computations
            [t,x,xr,time,timer,timer_noerr] = rmodel.getTrajectories(mu, inputidx);
            if this.UseOutput && ~isempty(rmodel.FullModel.System.C)
                x = rmodel.FullModel.System.C.computeOutput(t,x,mu);
                xr = rmodel.FullModel.System.C.computeOutput(t,xr,mu);
                est = rmodel.ErrorEstimator.OutputError;
            else
                est = rmodel.ErrorEstimator.StateError;
            end
            e = sqrt(sum((x - xr).^2,1));
            xnorm = sqrt(sum(x.^2,1));
            erel = e./xnorm;
            estrel = est./xnorm;
            xrnorm = sqrt(sum(xr.^2,1));
            erelr = e./xrnorm;
            estrelr = est./xrnorm;
            
            %% System plot
            xrmin = xr-repmat(est,size(xr,1),1); xrplus = xr+repmat(est,size(xr,1),1);
            ymax = max([max(x(:)) max(xr(:)) max(xrmin(:)) max(xrplus(:))]);
            ymin = min([min(x(:)) min(xr(:)) min(xrmin(:)) min(xrplus(:))]);
            h = figure;
            if this.SingleFigures
                subplot(1,2,1);
            else
                pos = get(0,'MonitorPosition');
                set(h,'OuterPosition',pos(1,:));
                subplot(2,3,1);
            end
            plot(t,general.Utils.preparePlainPlot(x));
            xlabel('t');
            if this.UseOutput
                title(sprintf('The full system''s output (m=%d,time=%.3f)',size(x,1),time));
            else
                title(sprintf('The full system (d=%d,time=%.3f)',size(x,1),time));
            end
            axis([t(1) t(end) ymin ymax]);
            if this.SingleFigures
                subplot(1,2,2);
            else
                subplot(2,3,2);
            end
            if this.UseOutput
                plot(t,general.Utils.preparePlainPlot(x),'b',...
                    t,general.Utils.preparePlainPlot(xr),'r',...
                    t,general.Utils.preparePlainPlot(xrmin),...
                    'r--',t,general.Utils.preparePlainPlot(xrplus),'r--');
                title(sprintf('Full+reduced system outputs with error bounds (r=%d,self time:%.3f, time with err est:%.3f)',size(rmodel.V,2),timer_noerr,timer));
                ylabel('y(t) / y^r(t)');
                legend('Full system','Reduced system','Lower bound','Upper bound');
            else
                plot(t,general.Utils.preparePlainPlot(xr),'r');
                title(sprintf('Reduced system (r=%d,self time:%.3f, time with err est:%.3f)',size(rmodel.V,2),timer_noerr,timer));
                ylabel('x^r(t)');
            end
            xlabel('t'); 
            axis([t(1) t(end) ymin ymax]);
            
            
            %% Absolute value plot
            if this.SingleFigures
                figure;
                subplot(1,2,1);
            else
                subplot(2,3,3);
            end
            plot(t,xnorm,'b',t,xrnorm,'r',t,xrnorm-est,'r--',t,xrnorm+est,'r--');
            xlabel('t');
            title('The state variable norms');
            legend('Full system','Reduced system','Lower bound','Upper bound');
            
            % Error plots
            if this.SingleFigures
                subplot(1,2,2);
            else
                subplot(2,3,4);
            end
            semilogy(t,est,'b',t,e,'r');%,t,abs(e-est),'g');
            xlabel('t');
            title(sprintf('The state variable absolute errors.\nmean(e)=%g, mean(est)=%g',mean(e),mean(est)));
            legend('Estimated error','True error');%,'Location','Best');
            
            % Relative Error plots
            if this.SingleFigures
                figure;
                subplot(1,2,1);
            else
                subplot(2,3,5);
            end
            semilogy(t,estrel,'b',t,erel,'r');
            xlabel('t');
            title(sprintf(['The state variable relative errors (comp. to '...
                'full solution)\nmean(e_{rel})=%g, mean(est_{rel})=%g'],...
                mean(erel),mean(estrel)));
            legend('Estimated error','True error');%,'Location','Best');
            
            if this.SingleFigures
                subplot(1,2,2);
            else
                subplot(2,3,6);
            end            
            semilogy(t,erelr,'r',t,estrelr,'b');
            xlabel('t');
            title(sprintf(['The state variable relative errors (comp. to '...
                'reduced solution)\nmean(ered_{rel})=%g, mean(estred_{rel})=%g'],...
                mean(erelr),mean(estrelr)));
            legend('True error','Estimated error');%,'Location','Best');
        end
    end
end

