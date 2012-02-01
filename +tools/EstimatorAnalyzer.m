classdef EstimatorAnalyzer < handle
% Analysis class for the error estimators.
%
% The constructor either takes an existing model and uses this for the
% error estimator analysis or the input argument is a dimension number for
% the standard model used.
%
% @author Daniel Wirtz @date 2010-08-01
%
% @change{0,6,dw,2011-12-05} Moved this class from the \c demos/EstimatorDemo to tools.EstimatorAnalyzer 
%
% @change{0,4,dw,2011-05-20} Adopted to the new strategy pattern implemented for the
% LocalLipschitzFcn inside the LocalLipschitzErrorEstimator (now having a class instead of a function handle).
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The used model
        %
        % @type models.BaseFullModel
        Model;
        
        % The reduced model
        %
        % @type models.ReducedModel
        ReducedModel;
        
        % How many estimator iterations should be performed?
        %
        % Is an integer row vector.
        %
        % @type rowvec
        EstimatorIterations = [1 2];
        
        % Chooses the estimator versions. Set to 1 for use, 0 for not use.
        % 1: True Error
        % 2: GLE, Global Lipschitz estimator.
        % 3: LGL, Local Gradient Lipschitz
        % 4: LSL, Local Secant Lipschitz
        % 5: ILSL, Improved Local Secant Lipschitz
        % 6: LGL TD, Time Discrete Local Gradient Lipschitz
        % 7: LSL TD, Time Discrete Local Secant Lipschitz
        % 8: ILSL TD, Time Discrete Improved Local Secant Lipschitz
        % 9: Expensive best-estimator (with full traj simulation)
        %
        % @type rowvec
        EstimatorVersions = [1 1 0 0 1 0 0 1 1];
        
        % Whether to put the estimations and computations times into a
        % single figure using subplots or using separate figures.
        %
        % @type logical
        SingleFigures = false;
        
        % Determines whether to plot the errors on a logarithmic scale or
        % not
        %
        % @type logical
        LogarithmicPlot = true;
        
        % Flag whether to save the time results as a LaTeX-Table into a
        % text file
        %
        % Set to [] to disable
        %
        % @type char
        SaveTexTables = 'tables.tex';
        
        % Set flag to sort the resulting computation time and estimates
        % table by the product of comp-time and error estimate.
        %
        % @type logical
        SortResultTable = false;
        
        % Flag wether to use the output errors or state variable errors
        %
        % @type logical
        UseOutputError = true;
        
        % The figure handles for the created figures.
        %
        % @type handle
        Figures;
        
        % Axes
        %
        % @type handle
        Axes;
        
        % The marker size for error and relative errors plots.
        MarkerSize = 8;
        
        % The number of markers for error and relative errors plots.
        NumMarkers = 5;
    end
    
    properties(SetAccess=private)
        % A struct containing information about different error estimators.
        %
        % @type struct
        Est;

        % A struct with fields Name, ErrsT, RelErrsT that contains
        % information about the estimations of the different error
        % estimators.
        ModelData;
    end
    
    methods
        
        function this = EstimatorAnalyzer(model)
            this.ModelData = struct('Name',{},'ErrsT',{},'RelErrsT',{});
            if nargin == 1
                this.setModel(model);
            end
        end
        
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
        
        function [errs, ctimes, varargout] = getErrorEstimates(this, mu, inidx, withrel)
            if nargin == 3
                withrel = false;
            end
            num = length(this.Est);
            ctimes = zeros(1,num);
            errs = zeros(num,length(this.Model.Times));
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulations
            %fprintf('Starting estimator demo for model %s!\n',this.Model.Name);
            % Save old estimator
            oldest = this.ReducedModel.ErrorEstimator;
            str = ''; compplot = [];
            for idx = 1:num
                % Set estimator and run
                this.ReducedModel.ErrorEstimator = this.Est(idx).Estimator;
                
                fprintf('Performing estimation %d of %d: %s...\n',idx,num,this.Est(idx).Name);
                [~, ~, ctimes(idx)] = this.ReducedModel.simulate(mu,inidx);
                if this.UseOutputError
                    errs(idx,:) = this.ReducedModel.ErrorEstimator.OutputError;
                else
                    errs(idx,:) = this.ReducedModel.ErrorEstimator.StateError;
                end
                % Plotting preparations
                if ~isa(this.Est(idx).Estimator,'error.ExpensiveBetaEstimator')
                    str = [str sprintf('errs(%d,end),ctimes(%d),''%s'',',idx,idx,this.Est(idx).MarkerStyle)]; %#ok<*AGROW>
                    compplot(end+1) = idx;
                end
            end
            % Restore old estimator
            this.ReducedModel.ErrorEstimator = oldest;
            
            if withrel
                varargout{1} = this.getRelativeErrorEstimates(errs, mu, inidx);
            end
        end
        
        function relerrs = getRelativeErrorEstimates(this, errs, mu, inidx)
            if this.UseOutputError
                [~,yf] = this.Model.simulate(mu, inidx);
            else
                [~,yf] = this.Model.computeTrajectory(mu, inidx);
                yf = bsxfun(@times, yf, this.ReducedModel.System.StateScaling);
            end
            yfullnorm = sqrt(sum(yf.*(this.Model.G*yf),1));
            relerrs = errs ./ repmat(yfullnorm,size(errs,1),1);
        end
        
        function plotErrors(this, errs)
            this.Figures{1} = figure;
            this.Axes{1} = gca;
            if ~this.SingleFigures
                pos = get(0,'MonitorPosition');
                set(this.Figures{1},'OuterPosition',pos(1,:));
                subplot(1,3,1);
            end
            if this.LogarithmicPlot
                ph = semilogy(this.Model.Times,errs);
            else
                ph = plot(this.Model.Times,errs);
            end
            set(ph(1),'LineWidth',2);
            hold on;
            % Select extra marker places
            nt = length(this.Model.Times);
            sel = round(1:nt/this.NumMarkers:nt);
            % Add some markers 
            for idx=1:length(this.Est)
                set(ph(idx),'LineStyle',this.Est(idx).LineStyle);
                % Shift marker positions for better visual
                pos = mod(sel+round(nt/(this.NumMarkers*length(this.Est)))*(idx-1),nt)+1;
                h = plot(this.Model.Times(pos), errs(idx,pos),this.Est(idx).MarkerStyle,'MarkerSize',this.MarkerSize);
                c = get(ph(idx),'Color');
                set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c*.7);
            end
            %keyboard;
            xlabel('Time');
            ylabel('Error estimates');
            a = cell(1,length(this.Est));
            [a{:}] = this.Est(:).Name;
            [~,oh] = legend(a,'Location','NorthEast');
            % Assign markers to legend
            oh = findobj(oh,'Type','line');
            for idx=1:length(this.Est)
                % Why ever, there are 2*numPlots line handles, and the second one controls the
                % markers in the center..
                c = get(oh(2*idx),'Color');
                set(oh(2*idx),'Marker',this.Est(idx).MarkerStyle,'MarkerFaceColor',c,'MarkerEdgeColor',c*.7,...
                    'MarkerSize',this.MarkerSize);
            end
            title(['Error estimations for model: ' this.Model.Name]);
            emin = min(errs(:));
            axis([0 this.Model.T emin*.9 max(emin,1e4)]);
        end
        
        function plotRelativeErrors(this, relerrs)
            if ~this.SingleFigures
                if isempty(this.Figures)
                    this.Figures{1} = figure;
                    pos = get(0,'MonitorPosition');
                    set(this.Figures{1},'OuterPosition',pos(1,:));
                end
                subplot(1,3,2);
            else
                this.Figures{2} = figure;
                this.Axes{2} = gca;
            end
            if this.LogarithmicPlot
                ph = semilogy(this.Model.Times,relerrs);
            else
                ph = plot(this.Model.Times,relerrs);
            end
            set(ph(1),'LineWidth',2);
            hold on;
            % Select extra marker places
            nt = length(this.Model.Times);
            sel = round(1:nt/this.NumMarkers:nt);
            for idx=1:length(this.Est)
                set(ph(idx),'LineStyle',this.Est(idx).LineStyle);
                % Shift marker positions for better visual
                pos = mod(sel+round(nt/(this.NumMarkers*length(this.Est)))*(idx-1),nt)+1;
                h = plot(this.Model.Times(pos), relerrs(idx,pos),this.Est(idx).MarkerStyle,...
                    'MarkerSize',this.MarkerSize);
                c = get(ph(idx),'Color');
                set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c*.7);
            end
            xlabel('Time');
            ylabel('Relative error estimates');
            a = cell(1,length(this.Est));
            [a{:}] = this.Est(:).Name;
            legend(a,'Location','NorthWest');
            [~,oh] = legend(a,'Location','NorthEast');
            % Assign markers to legend
            oh = findobj(oh,'Type','line');
            for idx=1:length(this.Est)
                % Why ever, there are 2*numPlots line handles, and the second one controls the
                % markers in the center..
                c = get(oh(2*idx),'Color');
                set(oh(2*idx),'Marker',this.Est(idx).MarkerStyle,'MarkerFaceColor',c,...
                    'MarkerEdgeColor',c*.7,'MarkerSize',this.MarkerSize);
            end
            title(['Relative error estimations e(t)/||y||, \Delta x(t)/||y|| for model: ' this.Model.Name]);
            emin = min(relerrs(:));
            axis([0 this.Model.T emin*.9 max(emin,1)]);
        end
        
        function plotCTimes(this, errs, ctimes)
            if ~this.SingleFigures
                if isempty(this.Figures)
                    this.Figures{1} = figure;
                    pos = get(0,'MonitorPosition');
                    set(this.Figures{1},'OuterPosition',pos(1,:));
                end
                subplot(1,3,3);
            else
                this.Figures{3} = figure;
                this.Axes{3} = gca;
            end
            str = ''; compplot = [];
            for idx = 1:length(this.Est)
                % Plotting preparations
                if ~isa(this.Est(idx).Estimator,'error.ExpensiveBetaEstimator')
                    str = [str sprintf('errs(%d,end),ctimes(%d),''%s'',',idx,idx,this.Est(idx).MarkerStyle)]; %#ok<*AGROW>
                    compplot(end+1) = idx;
                end
            end
            if this.LogarithmicPlot
                eval(['ph = semilogx(' str '''MarkerSize'',10);']);
            else
                eval(['ph = plot(' str '''MarkerSize'',10);']);
            end
            % Fill symbols
            for i=1:length(ph) 
                set(ph(i),'MarkerFaceColor',get(ph(i),'Color')); 
                set(ph(i),'MarkerEdgeColor',get(ph(i),'Color')*.7); 
            end
            % Add line for minimum error!
            hold on;
            plot([errs(1,end) errs(1,end)],[min(ctimes(compplot)) max(ctimes(compplot))],'black');
            hold off;
            
            xlabel('\Delta(T)');
            ylabel('Comp. time [s]');
            a = cell(1,length(this.Est));
            [a{:}] = this.Est(:).Name;
            legend(a(compplot));
            title(['Error estimator computation times: ' this.Model.Name]);
            axis tight;
        end
        
        function pt = getResultTable(this, errs, ctimes)
            if this.SortResultTable
                str = 'Estimator hierarchy (error*time product)';
                [~,idx] = sort(errs(:,end).*ctimes');
            else
                str = 'Estimator list';
                idx = 1:size(errs,2);
            end
            pt = PrintTable;
            pt.Caption = sprintf('%s for model "%s"',str,this.Model.Name);
            pt.addRow('Name',sprintf('$\\Delta(%d)$',this.Model.T),'Time','Overestimation');
            pt.HasHeader = true;
            for id = 1:length(this.Est)
                pt.addRow(this.Est(idx(id)).Name,errs(idx(id),end),ctimes(idx(id)),...
                    errs(idx(id),end)/errs(1,end),{'%s','$%1.3e}$','%2.2fs','$%1.3e}$'});
            end
        end
        
        function [ctimes, errs, relerrs] = start(this, mu, inidx)
            % Runs the demo with the current settings.
            %
            % Parameters:
            % mu: The parameter `\mu` to use @type colvec
            % inidx: The input index `i` of the input function `u_i` to use
            % @type integer
            
            if nargin < 3
                inidx = [];
                if nargin < 2
                    mu = [];
                end
            end
            
            [errs, ctimes, relerrs] = this.getErrorEstimates(mu, inidx, true);
            
            % Plot preps
            fprintf('Preparing plots...\n');
            this.Figures = {};
            
            % Absolute error plots
            this.plotErrors(errs)
            
            % Relative error plots
            this.plotRelativeErrors(relerrs);
            
            % Computation ctimes plot
            this.plotCTimes(errs, ctimes);
            
            % MaxErr data
            this.ModelData(end+1).Name = this.Model.Name;
            this.ModelData(end).ErrT = errs(:,end)';
            this.ModelData(end).MinErr = min(errs(:));
            this.ModelData(end).OverestT = (errs(:,end)/errs(1,end))';
            this.ModelData(end).RelErrT = relerrs(:,end)';
            this.ModelData(end).MinRelErr = min(relerrs(:));
            this.ModelData(end).CTimes = ctimes;
            this.ModelData(end).mu = mu;
            this.ModelData(end).inputidx = inidx;
            
            % Table overview
            t = this.getResultTable(errs, ctimes);
            t.display;
            
            if ~isempty(this.SaveTexTables)
                t.Format = 'tex';
                t.saveToFile(this.SaveTexTables);
%                 fid = fopen(this.SaveTexTables,'a+');
%                 str = [strrep(strrep(strrep(pt.print,'e+','\\cdot10^{'),'e-',...
%                         '\\cdot10^{-'),'{0','{') '\\\\\n'];
%                 fprintf(fid,str);
%                 fclose(fid);
            end
        end
        
        function ts = createStatsTables(this, sort)
            % Creates LaTeX tables with each the 'Errors','Relative errors','Overestimations' and
            % 'Computation times' for the demos started since creation of the tools.EstimatorAnalyzer.
            %
            % Uses the current model's name and param/input values as identification.
            %
            % Appends the resulting LaTeX table to the file 'statsTables.txt', which is created if
            % not existent.
            m = length(this.ModelData);
            if nargin < 2
                sort = 1:m;
            end
            fields = {'ErrT','RelErrT','OverestT','CTimes'};
            fieldnames = {'Errors','Relative errors','Overestimations','Computation times'};
            
            for fi = 1:length(fields)
                t = PrintTable;
                t.Format = 'tex';
                t.HasHeader = true;
                t.Caption = sprintf('%s of estimation runs',fieldnames{fi});
                t.addRow('Model / Est', this.Est(:).Name);
                
                fmt = '$%1.1e$';
                if strcmp(fields{fi},'CTimes')
                    fmt = '$%2.2fs$';
                end
                fmt = [{'$%s$'} repmat({fmt},1,length(this.Est))];
                for it=1:m
                    md = this.ModelData(sort(it));
                    
                    %% Compile name and add
                    nstr = '';
                    if ~isempty(md.mu)
                        nstr = [', \\mu=[' sprintf('%2.2e ',md.mu) ']'];
                    end
                    if ~isempty(md.inputidx)
                        nstr = [nstr sprintf(', u_%d',md.inputidx)];
                    end
                    hlp = num2cell(md.(fields{fi}));
                    t.addRow(nstr,hlp{:},fmt);
                end
                %str = strrep(strrep(strrep(strrep(t.print,'e+','e^{'),'e-','e^{-'),'{0','{'),'{-0','{-');
                %fprintf(fid,'%s',str);
                %t.print(fid);
                ts(fi) = t;
            end
            %fclose(fid);
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
                k = this.Model.Approx.Kernel;
            else
                k = this.Model.System.f.Kernel;
            end
            est = struct.empty;
            
            if this.EstimatorVersions(1)
                est(end+1).Name = 'True error';
                est(end).Estimator = error.DefaultEstimator(r);
                est(end).Estimator.Enabled = true;
                est(end).MarkerStyle = 'o';
                est(end).LineStyle = '-';
            end
            
            if this.EstimatorVersions(2)
                msg = error.GLEstimator.validModelForEstimator(r);
                if isempty(msg)
                        fprintf('Initializing Global Lipschitz estimator...\n');
                        est(end+1).Name = 'GLE';
                        est(end).Estimator = error.GLEstimator(r);
                        est(end).MarkerStyle = 's';
                        est(end).LineStyle = '-';
                else
                    fprintf('Cannot use the GLEstimator for model %s:\n%s\n',r.Name,msg);
                end
            end
            
            msg = error.IterationCompLemmaEstimator.validModelForEstimator(r);
            if isempty(msg)
                reps = this.EstimatorIterations;
                fprintf('Using iteration counts: %s\n',num2str(this.EstimatorIterations));
                
                if this.EstimatorVersions(3)
                    fprintf('Initializing LGL estimator...\n');
                    est(end+1).Name = 'LGL';
                    est(end).Estimator = error.IterationCompLemmaEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = error.lipfun.LocalGradientLipschitz(k);
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
                
                if this.EstimatorVersions(4)
                    fprintf('Initializing LGLMod (mod secant) estimator...\n');
                    est(end+1).Name = 'LGLMod';
                    est(end).Estimator = error.IterationCompLemmaEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = error.lipfun.LocalSecantLipschitz(k);
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
                
                if any(this.EstimatorVersions([5 8]))
                    ilfcn = error.lipfun.ImprovedLocalSecantLipschitz(k);
                end
                
                if this.EstimatorVersions(5)
                    fprintf('Initializing LSL estimator...\n');
                    est(end+1).Name = 'LSLE';
                    est(end).Estimator = error.IterationCompLemmaEstimator(r);
                    est(end).Estimator.LocalLipschitzFcn = ilfcn.clone;
                    est(end).Estimator.UseTimeDiscreteC = false;
                    est(end).MarkerStyle = 'p';
                    est(end).LineStyle = '-';
                    for it = reps
                        orig = est(end).Estimator;%#ok
                        eval(sprintf(['est(end+1).Name = ''LSLE, %d It'';'...
                            'est(end).Estimator = orig.clone;'...
                            'est(end).Estimator.Iterations = %d;'...
                            'est(end).MarkerStyle = ''p'';'...
                            'est(end).LineStyle = ''--'';'],it,it));
                    end
                end
                
                %% TD-Versions
                if any(this.EstimatorVersions(6:8))
                    td = error.IterationCompLemmaEstimator(r);
                    td.Iterations = 0;
                    td.UseTimeDiscreteC = true;
                end
                if this.EstimatorVersions(6)
                    fprintf('Initializing LSL TD estimators...\n');
                    est(end+1).Name = 'LGL TD';
                    est(end).Estimator = td;
                    est(end).Estimator.LocalLipschitzFcn = error.lipfun.LocalGradientLipschitz(k);
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                end
                
                if this.EstimatorVersions(7)
                    fprintf('Initializing LSL TD estimators...\n');
                    est(end+1).Name = 'LGL TD';
                    est(end).Estimator = td.clone;
                    est(end).Estimator.LocalLipschitzFcn = error.lipfun.LocalSecantLipschitz(k);
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                end
                
                if this.EstimatorVersions(8)
                    est(end+1).Estimator = est(end).Estimator.clone;
                    est(end).Name = 'LSLE TD';
                    est(end).Estimator = td.clone;
                    est(end).Estimator.LocalLipschitzFcn = ilfcn.clone;
                    est(end).MarkerStyle = '<';
                    est(end).LineStyle = '-';
                end
            else
                fprintf('Cannot use the IterationCompLemmaEstimator for model %s:\n%s\n',r.Name,msg);
            end
            
            if this.EstimatorVersions(9)
                msg = error.ExpensiveBetaEstimator.validModelForEstimator(r);
                if isempty(msg)
                        fprintf('Initializing expensive estimators with custom beta ...\n');
                        est(end+1).Name = 'Lower bound';
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

