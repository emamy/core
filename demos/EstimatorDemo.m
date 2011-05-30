classdef EstimatorDemo < handle
    % Demo class for the error estimators.
    %
    % The constructor either takes an existing model and uses this for the
    % error estimator demo or the input argument is a dimension number for
    % the standard model used.
    %
    % @change{0,4,dw,2011-05-20} Adopted to the new strategy pattern implemented for the
    % LocalLipschitzFcn inside the LocalLipschitzErrorEstimator (now having a class instead of a function handle).
    
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
        % 6: ILSL TD, experimental (with full traj simulation)
        EstimatorVersions = [1 0 0 0 1 1];
        
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
        
        % Flag wether to use the output errors or state variable errors
        UseOutputError = true;
    end
    
    properties(Access=private)
        % estimator struct
        Est;
    end
    
    methods
        
        function setModel(this, model)
            % Sets the model to use for the estimator demo
            if isa(model,'models.ReducedModel')
                this.Model = model.FullModel;
                this.ReducedModel = model;
            else
                this.Model = model;
            
                model.offlineGenerations;
                this.ReducedModel = model.buildReducedModel;
            end
            this.buildEstimatorStruct(this.ReducedModel);
        end
        
        function [ctimes, errs] = start(this, mu, inidx)
            % Runs the demo with the current settings.
            %
            % Texty text.
            
            if nargin < 3
                inidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            
            num = length(this.Est);
            ctimes = zeros(num,1);
            nt = length(this.Model.Times);
            errs = zeros(num,nt);
            
            % Plotting preparations
            str = '';
            h = waitbar(0,'');
            for idx = 1:num
                % Set estimator and run
                this.ReducedModel.ErrorEstimator = this.Est(idx).Estimator;
                
                waitbar(idx/num,h,sprintf('Performing estimation %d of %d: %s',idx,num,this.Est(idx).Name));
                [t,y,ctimes(idx)] = this.ReducedModel.simulate(mu,inidx);
                if this.UseOutputError
                    errs(idx,:) = this.ReducedModel.ErrorEstimator.OutputError;
                else
                    errs(idx,:) = this.ReducedModel.ErrorEstimator.StateError;
                end
                
                % Plotting preparations
                str = [str sprintf('errs(%d,end),ctimes(%d),''%s'',',idx,idx,this.Est(idx).MarkerStyle)];%#ok
            end
            close(h);
            % Compute full solution & remember full system's trajectory
            % (for rel. errors)
            %xrfullnorm = sqrt(sum(xr.^2,1));
            if this.UseOutputError
                [t,y] = this.Model.simulate(mu, inidx);
            else
                [t,y] = this.Model.computeTrajectory(mu, inidx);
                y = bsxfun(@times, y, this.ReducedModel.System.StateScaling);
                t = t*this.ReducedModel.tau;
            end
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
            sel = round(linspace(1,nt,6));
            sel = sel(1:5);
            
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
                pos = mod(sel+2*(idx-1),nt)+1;
                h = plot(this.Model.Times(pos), errs(idx,pos),this.Est(idx).MarkerStyle);
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
            relerrs = errs ./ repmat(yfullnorm,num,1);
            if this.LogarithmicPlot
                ph = semilogy(this.Model.Times,relerrs);
            else
                ph = plot(this.Model.Times,relerrs);
            end
            set(ph(1),'LineWidth',2);
            hold on;
            for idx=1:length(this.Est)
                set(ph(idx),'LineStyle',this.Est(idx).LineStyle);
                pos = mod(sel+2*(idx-1),nt)+1;
                h = plot(this.Model.Times(pos), relerrs(idx,pos),this.Est(idx).MarkerStyle);
                c = get(ph(idx),'Color');
                set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c*.7);
            end
            xlabel('Time');
            ylabel('Relative error estimates');
            legend(a,'Location','NorthWest');
            title(['Relative error estimations e(t)/||y||, \Delta x(t)/||y|| for model: ' this.Model.Name]);
            axis([0 this.Model.T 0 2]);
            
            %% Computation ctimes plot
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
            plot([errs(1,end) errs(1,end)],[min(ctimes) max(ctimes)],'black');
            hold off;
            
            xlabel('\Delta(T)');
            ylabel('Comp. time [s]');
            legend(a);
            title(['Error estimator computation times: ' this.Model.Name]);
            try
            axis tight;
            catch ME
            end
            
            %% Table overview
            if this.SortResultTable
                disp('Estimator hierarchy (error*time product):');
                [v,idx] = sort(errs(:,end).*ctimes);
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
                    errs(idx(id),end),ctimes(idx(id)),errs(idx(id),...
                    end)/errs(1,end),...
                    this.Est(idx(id)).Name);
                if ~isempty(this.SaveTexTables)
                    str = sprintf('\t\t%s & $%1.3e}$ & %2.2fs & $%1.3e}$',...
                        this.Est(idx(id)).Name,...
                        errs(idx(id),end),ctimes(idx(id)),...
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
                this.buildEstimatorStruct(this.ReducedModel);%#ok
            end
        end
        
        function set.EstimatorVersions(this, value)
            this.EstimatorVersions = value;
            if ~isempty(this.ReducedModel)%#ok
                this.buildEstimatorStruct(this.ReducedModel);%#ok
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
            
            if this.EstimatorVersions(1)
                msg = error.GLEstimator.validModelForEstimator(r);
                if isempty(msg)
                        fprintf('Initializing Global Lipschitz estimator...\n');
                        est(end+1).Name = 'GLE';
                        est(end).Estimator = error.GLEstimator(r);
                        est(end).MarkerStyle = 'o';
                        est(end).LineStyle = '-';
                else
                    fprintf('Cannot use the GLEstimator for model %s:\n%s\n',r.Name,msg);
                end
            end
            
            msg = error.LocalKernelEstimator.validModelForEstimator(r);
            if isempty(msg)
                reps = this.EstimatorIterations;
                fprintf('Using iteration counts: %s\n',num2str(this.EstimatorIterations));
                
                if this.EstimatorVersions(2)
                    fprintf('Initializing LGL estimator...\n');
                    est(end+1).Name = 'LGL';
                    est(end).Estimator = error.LocalKernelEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = error.LocalGradientLipschitz(k);
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
                    est(end).Estimator = error.LocalKernelEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = error.LocalSecantLipschitz(k);
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
                
                if any(this.EstimatorVersions(4:5))
                    ilfcn = error.ImprovedLocalSecantLipschitz(k);
                end
                
                if this.EstimatorVersions(4)
                    fprintf('Initializing LSL estimator...\n');
                    est(end+1).Name = 'LSL';
                    est(end).Estimator = error.LocalKernelEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = ilfcn;
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
                    est(end).Estimator = error.LocalKernelEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = error.LocalGradientLipschitz(k);
                    est(end).Estimator.UseTimeDiscreteC = true;
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                    
                    est(end+1).Estimator = est(end).Estimator.clone;
                    est(end).Name = 'LSL TD';
                    est(end).Estimator.LocalLipschitzFcn = ilfcn;
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                end
            else
                fprintf('Cannot use the LocalKernelEstimator for model %s:\n%s\n',r.Name,msg);
            end
            
            if this.EstimatorVersions(6)
                msg = error.ExpensiveBetaEstimator.validModelForEstimator(r);
                if isempty(msg)
                        fprintf('Initializing expensive estimators with custom beta ...\n');
                        est(end+1).Name = 'TD Bestest';
                        est(end).Estimator = error.ExpensiveBetaEstimator(r);
                        est(end).MarkerStyle = '<';
                        est(end).LineStyle = '-.';

%                         est(end+1).Estimator = est(end).Estimator.clone;
%                         est(end).Name = 'TD Jacobian (nonrigor)';
%                         est(end).Estimator.Version = 2;
%                         est(end).MarkerStyle = '<';
%                         est(end).LineStyle = '-.';
                else
                    fprintf('Cannot use the ExpensiveBetaEstimator for model %s:\n%s\n',r.Name,msg);
                end
            end
            this.Est = est;
        end
    end
end

