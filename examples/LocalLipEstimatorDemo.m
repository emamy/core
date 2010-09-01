classdef LocalLipEstimatorDemo < handle
    % Demo class for the error estimators.
    %
    % The constructor either takes an existing model and uses this for the
    % error estimator demo or the input argument is a dimension number for
    % the standard model used.
    
    properties
        % System dimension
        Dims = 500;
        
        % Number of centers used in kernel expansion
        NumCenters = 10;
        
        % Strictly positive kernel expansion?
        PositiveExpansion = false;
        
        % Uniform expansion or comp-wise separate?
        UniformExpansion = false;
        
        % The used model
        Model;
        
        % Reduced model
        r;
    end
    
    properties(Access=private)
        % estimator struct
        Est;
    end
    
    methods
        
        function this = LocalLipEstimatorDemo(inarg)
            % Creates a new estimator demo.
            %
            % Parameters:
            % inarg: Is either a models.BaseFullModel subclass that
            % determines the model to use OR a scalar defining the
            % dimensions if the default test model. Calling the constructor
            % with no arguments causes the Demo to use the default demo
            % model with the default dimension size (500).
            
            % Catch the model argument case
            if nargin == 1 && isa(inarg,'models.BaseFullModel')
                this.setModel(inarg);
                return;
            end
            
            %% Model settings
            fm = models.BaseFullModel;
            fm.Verbose = 0;
            fm.Name = 'Estimator Demo Model';
            
            fm.T = 1;
            fm.dt = .025;
            
            fm.Approx = [];
            fm.Sampler = [];
            
            %this.ODESolver = solvers.MLWrapper(@ode45);
            fm.ODESolver = solvers.ExplEuler(fm.dt);
            %fm.ODESolver = solvers.Heun(fm.dt);
            
            %% Core function
            cf = dscomponents.CompwiseKernelCoreFun;
            cf.TimeKernel = kernels.NoKernel;
            cf.ParamKernel = kernels.NoKernel;
            cf.snData.ti = [];
            cf.snData.mui = [];
            
            %% System settings
            fm.System = models.BaseDynSystem;
            fm.System.f = cf;
            this.Model = fm;
            
            if nargin == 1
                this.Dims = inarg;
            else
                this.newCoeffs;
                %this.setup;
            end
        end
        
        function setup(this)
            this.Model.System.f.SystemKernel = kernels.GaussKernel(15);
            x0 = rand(this.Dims,1);
            if this.PositiveExpansion
                base = linspace(0, 40, this.NumCenters);
                this.Model.System.x0 = @(mu)x0;
            else
                base = linspace(-20, 20, this.NumCenters);
                this.Model.System.x0 = @(mu)x0-.5;
            end
            this.Model.System.f.snData.xi = repmat(base,this.Dims,1);
            
            if  this.UniformExpansion
                V = ones(this.Dims,1)*sqrt(1/this.Dims);
                s = spacereduction.ManualReduction(V,V);
            else
                s = spacereduction.PODReducer;
                s.Mode = 'abs';
                s.Value = 1;
            end
            this.Model.SpaceReducer = s;
            
            %% Generation
            this.Model.offlineGenerations;
            this.r = this.Model.buildReducedModel;
            
            this.buildEstimatorStruct(this.r);
        end
        
        function newCoeffs(this)
            % Function coefficients
            offset = .5;
            if this.PositiveExpansion
                offset = 0;
            end
            % Create coefficients
            if this.UniformExpansion
                ai = (rand(1,this.NumCenters)-offset);
                this.Model.System.f.Ma = repmat(ai,this.Dims,1);
            else
                this.Model.System.f.Ma = (rand(this.Dims,this.NumCenters)-offset);
            end
            this.setup;
        end
        
        function setModel(this, model)
            this.Model = model;
            
            model.offlineGenerations;
            this.r = model.buildReducedModel;
            
            this.buildEstimatorStruct(this.r);
        end
        
        function [times,errs] = Run(this, mu, inputidx)
            % Runs the demo with the current settings.
            if nargin < 3
                inputidx = 1;
                if nargin < 2
                    mu = [.2; .6; -.4];
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
                this.r.ErrorEstimator = this.Est(idx).Estimator;
                
                waitbar(idx/num,h,sprintf('Performing estimation %d of %d: %s',idx,num,this.Est(idx).Name));
                [t,y,times(idx)] = this.r.simulate(mu,inputidx);
                errs(idx,:) = this.r.ErrorEstimator.LastError;
                
                % Plotting preparations
                str = [str sprintf('errs(%d,end),times(%d),''%s'',',idx,idx,this.Est(idx).MarkerStyle)];%#ok
            end
            close(h);
            
            % Cut global estimation error
%             tmp = errs(2,:);
%             tmp(tmp > errs(3,end)) = errs(3,end);
%             errs(2,:) = tmp;
            
            %% Plot
            h = figure;
            pos = get(0,'MonitorPosition');
            set(h,'OuterPosition',pos(1,:));
            a = cell(1,num);
            [a{:}] = this.Est(:).Name;
            
            subplot(1,2,1);
            ph = plot(this.Model.Times,errs);
            set(ph(1),'LineWidth',2);
            for idx=1:length(this.Est)
                set(ph(idx),'LineStyle',this.Est(idx).LineStyle);
            end
            xlabel('Time');
            ylabel('Error');
            legend(a,'Location','NorthWest');
            
            subplot(1,2,2);
            eval(['ph = plot(' str '''MarkerSize'',7);']);
            %set(ph(6:9),'Marker','+');
            %set(ph(10:11),'Marker','*');
            
            % Add line for minimum error!
            hold on;
            plot([errs(1,end) errs(1,end)],[min(times) max(times)],'black');
            hold off;
            
            %txt = [repmat('(',num,1) num2str(times) repmat(', ',num,1) num2str(errs(:,end)) repmat('s)',num,1)];
            %text(times,errs(:,end),txt);
            xlabel('e(T)');
            ylabel('Comp. time');
            legend(a);
            
            title(['Error estimator demo for model: ' this.Model.Name]);
            disp('Computation times:');
            disp(times');
            disp('e(T):');
            disp(errs(:,end)');
            
            disp('Best five estimations:');
            [v,idx] = sort(errs(:,end));
            for id = 1:5
                fprintf('%s: e(T)=%f, %2.2fsec\n',this.Est(idx(id)).Name,errs(idx(id),end),times(idx(id)));
            end
        end
        
        function set.Dims(this, value)
            this.Dims = value;
            this.newCoeffs; %#ok
            %this.setup; %#ok
        end
        
        function set.NumCenters(this, value)
            this.NumCenters = value;
            this.newCoeffs; %#ok
            %this.setup; %#ok
        end
        
        function set.PositiveExpansion(this, value)
            this.PositiveExpansion = value;
            %this.setup; %#ok
            this.newCoeffs; %#ok
        end
        
        function set.UniformExpansion(this, value)
            this.UniformExpansion = value;
            this.newCoeffs; %#ok
            %this.setup; %#ok
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
            
%             msg = error.GlobalLipKernelEstimator.validModelForEstimator(r);
%             if isempty(msg)
%                 est(end+1).Name = 'GLE';
%                 est(end).Estimator = error.GlobalLipKernelEstimator(r);
%                 est(end).MarkerStyle = 'x';
%                 est(end).LineStyle = '-';
%             else
%                 fprintf('Cannot use the GlobalLipKernelEstimator for model %s:\n%s\n',r.Name,msg);
%             end
            
            msg = error.LocalLipKernelEstimator.validModelForEstimator(r);
            if isempty(msg)
                reps = [1 2];
                est(end+1).Name = 'LLE: LGL';
                est(end).Estimator = error.LocalLipKernelEstimator(r);
                est(end).Estimator.KernelLipschitzFcn = @k.getLocalGradientLipschitz;
                est(end).MarkerStyle = 's';
                est(end).LineStyle = '-';
                for it = reps
                    eval(sprintf(['est(end+1).Name = ''LLE: LGL, %d It'';'...
                        'est(end).Estimator = error.LocalLipKernelEstimator(r);'...
                        'est(end).Estimator.KernelLipschitzFcn = @k.getLocalGradientLipschitz;'...
                        'est(end).Estimator.Iterations = %d;'...
                        'est(end).MarkerStyle = ''s'';'...
                        'est(end).LineStyle = '':'';'],it,it));
                end
                
                est(end+1).Name = 'LLE: LSL';
                est(end).Estimator = error.LocalLipKernelEstimator(r);
                est(end).Estimator.KernelLipschitzFcn = @k.getLocalSecantLipschitz;
                est(end).MarkerStyle = '*';
                est(end).LineStyle = '-';
                for it = reps
                    eval(sprintf(['est(end+1).Name = ''LLE: LSL, %d It'';'...
                        'est(end).Estimator = error.LocalLipKernelEstimator(r);'...
                        'est(end).Estimator.KernelLipschitzFcn = @k.getLocalSecantLipschitz;'...
                        'est(end).Estimator.Iterations = %d;'...
                        'est(end).MarkerStyle = ''*'';'...
                        'est(end).LineStyle = ''-.'';'],it,it));
                end
                
                est(end+1).Name = 'LLE: ILSL';
                est(end).Estimator = error.LocalLipKernelEstimator(r);
                est(end).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;
                est(end).MarkerStyle = '+';
                est(end).LineStyle = '-';
                for it = reps
                    eval(sprintf(['est(end+1).Name = ''LLE: ILSL, %d It'';'...
                        'est(end).Estimator = error.LocalLipKernelEstimator(r);'...
                        'est(end).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;'...
                        'est(end).Estimator.Iterations = %d;'...
                        'est(end).MarkerStyle = ''+'';'...
                        'est(end).LineStyle = ''--'';'],it,it));
                end
                
                est(end+1).Name = 'LLE: ILSL TD';
                est(end).Estimator = error.LocalLipKernelEstimator(r);
                est(end).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;
                est(end).Estimator.UseTimeDiscreteC = true;
                est(end).MarkerStyle = '<';
                est(end).LineStyle = '-';
%                 for it = 1
%                     eval(sprintf(['est(end+1).Name = ''LLE: ILSL TD, %d It'';'...
%                         'est(end).Estimator = error.LocalLipKernelEstimator(r);'...
%                         'est(end).Estimator.KernelLipschitzFcn = @k.getImprovedLocalSecantLipschitz;'...
%                         'est(end).Estimator.UseTimeDiscreteC = true;'...
%                         'est(end).Estimator.Iterations = %d;'...
%                         'est(end).MarkerStyle = ''<'';'...
%                         'est(end).LineStyle = ''--'';'],it,it));
%                 end
            else
                fprintf('Cannot use the LocalLipKernelEstimator for model %s:\n%s\n',r.Name,msg);
            end
            
            this.Est = est;
        end
    end
end

