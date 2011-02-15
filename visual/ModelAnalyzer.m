classdef ModelAnalyzer < handle;
    %MODELANALYZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SingleFigures = false;
        
        UseOutput = true;
    end
    
    methods
        
        function this = ModelAnalyzer(rmodel)
            if nargin > 0
                this.analyze(rmodel);
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
            end
            e = sqrt(sum((x - xr).^2,1));
            est = rmodel.ErrorEstimator.LastError;
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
            plot(t,est,'b',t,e,'r');%,t,abs(e-est),'g');
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
            plot(t,estrel,'b',t,erel,'r');
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
            plot(t,erelr,'r',t,estrelr,'b');
            xlabel('t');
            title(sprintf(['The state variable relative errors (comp. to '...
                'reduced solution)\nmean(ered_{rel})=%g, mean(estred_{rel})=%g'],...
                mean(erelr),mean(estrelr)));
            legend('True error','Estimated error');%,'Location','Best');
        end
        
%         function compareCoreVsApprox(this)
%             num = 10;
%             d = rmodel.FullModel.Data;
%             
%             % Get data
%             sn = d.ApproxTrainData;
%             
%             factor = 4;
%             % Extend data
%             nz = factor*size(sn,2);
%             oldidx = 1:factor:nz;
%             
%             xi = sn(4:end,:);
%             XI(:,oldidx) = xi;
%             di = (xi(:,2:end)-xi(:,1:end-1))/factor;
%             di(:,end+1) = 0;
%             ti = sn(3,:);
%             TI(oldidx) = ti;
%             dti = (ti(2:end)-ti(1:end-1))/factor;
%             dti(end+1) = 0;
%             mui(oldidx) = sn(1,:);
%             %fx = zeros(size(d.ApproxfValues,1),size(d.ApproxfValues,2)*factor);
%             fx(:,oldidx) = d.ApproxfValues;
%             for fac = 1:factor-1
%                 newidx = oldidx+fac;
%                 % xi
%                 XI(:,newidx) = xi+di*(fac/(factor-1));
%                 % ti
%                 TI(newidx) = ti+dti*(fac/(factor-1));
%                 % mui
%                 mui(newidx) = sn(1,:);
%                 
%                 % Evaluate original function at middle points
%                 for idx=newidx
%                     fx(:,idx) = rmodel.FullModel.System.f.evaluate(XI(:,idx),TI(idx),d.getParams(mui(idx)));
%                 end
%             end
%             MUI = d.getParams(mui);
%             
%             afx = rmodel.FullModel.Approx.evaluate(XI,TI,MUI);
%             
%             err = sqrt(sum((fx-afx).^2));
%             total_err = sqrt(sum(err.^2))
%             
%             % Plotting
%             len = size(afx,2);
%             muvals = (mui / max(mui));
%             for idx = 1:num
%                 fxp = fx(idx,:);
%                 afxp = afx(idx,:);
%                 subplot(1,2,1);
%                 plot(1:len,fxp,'r',1:len,afxp,'b--',1:len,muvals*(max(fxp)-min(fxp)),'green',[oldidx; oldidx + eps],[0 0],'black+');
%                 subplot(1,2,2);
%                 plot(1:len,err,'r',1:len,muvals*(max(err)-min(err)),'green',[oldidx; oldidx + eps],[0 0],'black+');
%                 pause;
%             end
%             
%             %plot(1:len,fx,1:len,afx,'--');
%             
%             %             arfx = rmodel.System.f.evaluate(d.V'*xi,ti,mui);
%             %             rfx = d.V'*d.ApproxfValues(:,:);
%             %             for idx = 1:num
%             %                 plot(1:len,rfx(idx,:),1:len,arfx(idx,:));
%             %                 pause;
%             %             end
%         end
        
        %         function save(this, matfile)
        %             % Saves the reduced model to disk.
        %             %
        %             % Parameters:
        %             % matfile: The target file. If not specified, a file with the
        %             % name of the reduced model's variable name is used.
        %
        %             name = inputname(1);
        %             if nargin == 1
        %                 matfile = fullfile(cd,name);
        %             end
        %             % For the save process of the reduced model the full model's
        %             % Data (=ModelData) and Approx properties are not needed. This
        %             % is the fastest way to ensure that the reduced model can still
        %             % have access to all important features of the full model but
        %             % uses less disk space.
        %             m = rmodel.FullModel;
        %             d = m.Data;
        %             a = m.Approx;
        %
        %             m.Data = [];
        %             m.Approx = [];
        %
        %             eval([name ' = this;']);
        %             save(matfile,name);
        %
        %             m.Data = d;
        %             m.Approx = a;
        %         end
    end
    
end

