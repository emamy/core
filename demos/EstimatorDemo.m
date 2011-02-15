classdef EstimatorDemo < handle
    % Demo class for the error estimators.
    %
    % The constructor either takes an existing model and uses this for the
    % error estimator demo or the input argument is a dimension number for
    % the standard model used.
    
    properties
        % The used model
        Model;
        
        % The reduced model
        ReducedModel;
        
        % How many estimator iterations should be performed?
        EstimatorIterations = [1 2];
        
        % Chooses the estimator versions. Set to 1 for use, 0 for not use.
        % 1: GLE, Global Lipschitz estimator.
        % 2: LGL, Local Gradient Lipschitz
        % 3: LSL, Local Secant Lipschitz
        % 4: ILSL, Improved Local Secant Lipschitz
        % 5: ILSL TD, Time Discrete Improved Local Secant Lipschitz
        EstimatorVersions = [1 1 0 1 1];
        
        % Whether to put the estimations and computations times into a
        % single figure using subplots or using separate figures.
        SingleFigures = false;
        
        % Determines whether to plot the errors on a logarithmic scale or
        % not
        LogarithmicPlot = true;
        
        % Flag whether to save the time results as a LaTeX-Table into a
        % text file
        %
        % Set to [] to disable
        SaveTexTables = 'tables.tex';
        
        % Set flag to sort the resulting computation time and estimates
        % table by the product of comp-time and error estimate.
        SortResultTable = false;
    end
    
    properties(Access=private)
        % estimator struct
        Est;
    end
    
    methods
        
        function setModel(this, model)
            % Sets the model to use for the estimator demo
            this.Model = model;
            
            model.offlineGenerations;
            this.ReducedModel = model.buildReducedModel;
            
            this.buildEstimatorStruct(this.ReducedModel);
        end
        
        function [times,errs] = Run(this, mu, inidx)
            % Runs the demo with the current settings.
            
            if nargin < 3
                inidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            
            num = length(this.Est);
            times = zeros(num,1);
            errs = zeros(num,length(this.Model.Times));
            
            % Plotting preparations
            str = '';
            h = waitbar(0,'');
            for idx = 1:num
                % Set estimator and run
                this.ReducedModel.ErrorEstimator = this.Est(idx).Estimator;
                
                waitbar(idx/num,h,sprintf('Performing estimation %d of %d: %s',idx,num,this.Est(idx).Name));
                [t,y,times(idx)] = this.ReducedModel.simulate(mu,inidx);
                errs(idx,:) = this.ReducedModel.ErrorEstimator.LastError;
                
                % Plotting preparations
                str = [str sprintf('errs(%d,end),times(%d),''%s'',',idx,idx,this.Est(idx).MarkerStyle)];%#ok
            end
            close(h);
            % Compute full solution & remember full system's trajectory
            % (for rel. errors)
            %xrfullnorm = sqrt(sum(xr.^2,1));
            [t,y] = this.Model.simulate(mu,inidx);
            yfullnorm = sqrt(sum(y.^2,1));
            %deferrest = errs(1,:);
            %trueerr = sqrt(sum((x-this.ReducedModel.V*xr).^2,1));
            
            %% Plot preps
            h = figure;
            a = cell(1,num);
            [a{:}] = this.Est(:).Name;
            if ~this.SingleFigures
                pos = get(0,'MonitorPosition');
                set(h,'OuterPosition',pos(1,:));
                subplot(1,3,1);
            end
            sel = round(linspace(1,length(this.Model.Times),6));
            sel = sel(2:5);
            
            %% Absolute error plots
            if this.LogarithmicPlot
                ph = semilogy(this.Model.Times,errs);
            else
                ph = plot(this.Model.Times,errs);
            end
            set(ph(1),'LineWidth',2);
            hold on;
            % Add some markers
            for idx=1:length(this.Est)
                set(ph(idx),'LineStyle',this.Est(idx).LineStyle);
                h = plot(this.Model.Times(sel+2*idx), errs(idx,sel+2*idx),this.Est(idx).MarkerStyle);
                c = get(ph(idx),'Color');
                set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c*.7);
            end
            %keyboard;
            xlabel('Time');
            ylabel('Error estimates');
            legend(a,'Location','NorthWest');
            title(['Error estimations for model: ' this.Model.Name]);
            axis tight;
            
            %% Relative error plots
            if ~this.SingleFigures
                subplot(1,3,2);
            else
                figure;
            end
            ph = plot(this.Model.Times,errs ./ repmat(yfullnorm,num,1));
            set(ph(1),'LineWidth',2);
            hold on;
            for idx=1:length(this.Est)
                set(ph(idx),'LineStyle',this.Est(idx).LineStyle);
                h = plot(this.Model.Times(sel+2*idx), errs(idx,sel+2*idx),this.Est(idx).MarkerStyle);
                c = get(ph(idx),'Color');
                set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c*.7);
            end
            xlabel('Time');
            ylabel('Relative error estimates');
            legend(a,'Location','NorthWest');
            title(['Relative error estimations e(t)/||y||, \Delta x(t)/||y|| for model: ' this.Model.Name]);
            axis([0 this.Model.T 0 2]);
            
            %% Computation times plot
            if ~this.SingleFigures
                subplot(1,3,3);
            else
                figure;
            end
            if this.LogarithmicPlot
                eval(['ph = semilogx(' str '''MarkerSize'',8);']);
            else
                eval(['ph = plot(' str '''MarkerSize'',8);']);
            end
            % Fill symbols
            for i=1:length(ph) 
                set(ph(i),'MarkerFaceColor',get(ph(i),'Color')); 
                set(ph(i),'MarkerEdgeColor',get(ph(i),'Color')*.7); 
            end
            % Add line for minimum error!
            hold on;
            plot([errs(1,end) errs(1,end)],[min(times) max(times)],'black');
            hold off;
            
            xlabel('\Delta(T)');
            ylabel('Comp. time [s]');
            legend(a);
            title(['Error estimator computation times: ' this.Model.Name]);
            axis tight;
            
            %% Table overview
            if this.SortResultTable
                disp('Estimator hierarchy (error*time product):');
                [v,idx] = sort(errs(:,end).*times);
            else
                disp('Estimator list:');
                idx = 1:size(errs,2);
            end
            if ~isempty(this.SaveTexTables)
                fid = fopen(this.SaveTexTables,'a+');
                fprintf(fid,'\n%% %s: Table for model %s\n',datestr(clock),this.Model.Name);
                fprintf(fid,'\\begin{table}[htbp]\n\t\\centering\\scriptsize\n\t\\begin{tabular}{l|r|r|l}\n');
                fprintf(fid,'\tName & $\\Delta(%d)$ & Time & Overestimation\\\\\\hline\n',t(end));
            end
            for id = 1:length(this.Est)
                fprintf('Delta(%2.2f)=%e\t%2.4fsec\tfactor:%e\t%s\n',t(end),...
                    errs(idx(id),end),times(idx(id)),errs(idx(id),...
                    end)/errs(1,end),...
                    this.Est(idx(id)).Name);
                if ~isempty(this.SaveTexTables)
                    str = sprintf('\t\t%s & $%1.3e}$ & %2.2fs & $%1.3e}$',...
                        this.Est(idx(id)).Name,...
                        errs(idx(id),end),times(idx(id)),...
                        errs(idx(id),end)/errs(1,end));
                    fprintf(fid,[strrep(strrep(strrep(str,'e+','\\cdot10^{'),'e-',...
                        '\\cdot10^{-'),'{0','{') '\\\\\n']);
                end
            end
            if ~isempty(this.SaveTexTables)
                fprintf(fid, '\t\\end{tabular}\n');
                fprintf(fid,'\t\\caption{Estimator hierarchy for model "%s"}\n',this.Model.Name);
                fprintf(fid, '\\end{table}\n');
                fprintf(fid,'%% End of output for model %s\n',this.Model.Name);
                fclose(fid);
            end
        end
        
        function set.EstimatorIterations(this, value)
            this.EstimatorIterations = value;
            if ~isempty(this.ReducedModel)%#ok
                this.buildEstimatorStruct(this, this.ReducedModel);%#ok
            end
        end
        
    end
    
    methods(Access=private)
        function buildEstimatorStruct(this, r)
            % Error estimators
            if ~isempty(this.Model.Approx)
                k = this.Model.Approx.SystemKernel;
            else
                k = this.Model.System.f.SystemKernel;
            end
            est = struct.empty;
            
            est(end+1).Name = 'Full error';
            est(end).Estimator = error.DefaultEstimator(r);
            est(end).Estimator.Enabled = true;
            est(end).MarkerStyle = 'o';
            est(end).LineStyle = '-';
            
            msg = error.GlobalLipKernelEstimator.validModelForEstimator(r);
            if isempty(msg)
                if this.EstimatorVersions(1)
                    fprintf('Initializing Global Lipschitz estimator...\n');
                    est(end+1).Name = 'GLE';
                    est(end).Estimator = error.GlobalLipKernelEstimator(r);
                    est(end).MarkerStyle = 'o';
                    est(end).LineStyle = '-';
                end
            else
                fprintf('Cannot use the GlobalLipKernelEstimator for model %s:\n%s\n',r.Name,msg);
            end
            
            msg = error.LocalLipKernelEstimator.validModelForEstimator(r);
            if isempty(msg)
                reps = this.EstimatorIterations;
                fprintf('Using iteration counts: %s\n',num2str(this.EstimatorIterations));
                    
                if this.EstimatorVersions(2)
                    fprintf('Initializing LGL estimator...\n');
                    est(end+1).Name = 'LGL';
                    est(end).Estimator = error.LocalLipKernelEstimator(r);
                    est(end).Estimator.KernelLipschitzFcn = @k.getLocalGradientLipschitz;
                    est(end).Estimator.UseTimeDiscreteC = false;
                    est(end).MarkerStyle = 's';
                    est(end).LineStyle = '-';
                    for it = reps
                        orig = est(end).Estimator;%#ok
                        eval(sprintf(['est(end+1).Name = ''LGL, %d It'';'...
                            'est(end).Estimator = orig.clone;'...
                            'est(end).Estimator.Iterations = %d;'...
                            'est(end).MarkerStyle = ''s'';'...
                            'est(end).LineStyle = ''-.'';'],it,it));
                    end
                end
                
                if this.EstimatorVersions(3)
                    fprintf('Initializing LSL (mod secant) estimator...\n');
                    est(end+1).Name = 'LGLMod';
                    est(end).Estimator = error.LocalLipKernelEstimator(r);
                    est(end).Estimator.KernelLipschitzFcn = @k.getLocalSecantLipschitz;
                    est(end).Estimator.UseTimeDiscreteC = false;
                    est(end).MarkerStyle = 'h';
                    est(end).LineStyle = '-';
                    for it = reps
                        orig = est(end).Estimator;%#ok
                        eval(sprintf(['est(end+1).Name = ''LGLMod, %d It'';'...
                            'est(end).Estimator = orig.clone;'...
                            'est(end).Estimator.Iterations = %d;'...
                            'est(end).MarkerStyle = ''h'';'...
                            'est(end).LineStyle = '':'';'],it,it));
                    end
                end
                
                if this.EstimatorVersions(4)
                    fprintf('Initializing LSL estimator...\n');
                    est(end+1).Name = 'LSL';
                    est(end).Estimator = error.LocalLipKernelEstimator(r);
                    est(end).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;
                    est(end).Estimator.UseTimeDiscreteC = false;
                    est(end).MarkerStyle = 'p';
                    est(end).LineStyle = '-';
                    for it = reps
                        orig = est(end).Estimator;%#ok
                        eval(sprintf(['est(end+1).Name = ''LSL, %d It'';'...
                            'est(end).Estimator = orig.clone;'...
                            'est(end).Estimator.Iterations = %d;'...
                            'est(end).MarkerStyle = ''p'';'...
                            'est(end).LineStyle = ''--'';'],it,it));
                    end
                end
                
                if this.EstimatorVersions(5)
                    fprintf('Initializing LSL TD estimators...\n');
                    est(end+1).Name = 'LGL TD';
                    est(end).Estimator = error.LocalLipKernelEstimator(r);
                    est(end).Estimator.KernelLipschitzFcn = @k.getLocalGradientLipschitz;
                    est(end).Estimator.UseTimeDiscreteC = true;
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                    
                    est(end+1).Estimator = est(end).Estimator.clone;
                    est(end).Name = 'LSL TD';
                    est(end).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                end
            else
                fprintf('Cannot use the LocalLipKernelEstimator for model %s:\n%s\n',r.Name,msg);
            end
            
            this.Est = est;
        end
    end
end

