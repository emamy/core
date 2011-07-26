classdef WH10Experiment < models.BaseFullModel
    % Numerical experiments class for Paper WH10
    %
    % Current version works with KerMor 0.4
    
    properties(SetObservable)
        % The system's dimension
        %
        % @propclass{experimental}
        dim;
    end
    
    methods
        
        function this = WH10Experiment(dims)
            this.registerProps('dim');
            
            this.dim = dims;
            
            this.Sampler = [];%sampling.GridSampler;
            
            % This class implements a fake Approx subclass to allow access
            % to the this.Ma property for the error estimator.
            this.Approx = [];
            
            %s = solvers.ode.MLWrapper(@ode45);
            s = solvers.ode.ExplEuler;
            s.MaxStep = [];
            %s = solvers.ode.Heun;
            this.ODESolver = s;
            
            %% System settings
            this.System = WH10System(this);
        end
    end    
    
    methods(Static)
        
        function [d, rmodels] = Experiment1(dim)
            % Experiment with manual reduction and rotation as error source
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
            f = m.System.f;
            
            f.x0 = dscomponents.ConstInitialValue(ones(dim,1));
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dim,1);
            
            %m.System.Inputs{1} = @(t).4*exp(-t/4); % alter input
            m.System.Inputs{1} = @(t).04*sin(t/3);
            m.System.Inputs{2} = @(t)exp(-abs(10-t/2)); % Gutes Beispiel!
            m.System.Inputs{3} = @(t).5*exp(-(12-t).^2); % Sehr interessante dynamik
            m.System.Inputs{4} = @(t)(t>10)*.1; % okay, aber nicht so spannend
            m.System.B = dscomponents.LinearInputConv(ones(m.dim,1));
            
            % No training necessary!
            m.TrainingInputs = [];
            
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            V = ones(dim,1)/sqrt(dim);
            s = spacereduction.ManualReduction(V,V);
            s = spacereduction.RotationDecorator(s);
            s.Dims = 100;
            m.SpaceReducer = s;
            
            d = EstimatorDemo;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            
            degs = [.0005 .005 .05 .5];
            rmodels = {};
            
            for didx = 1:length(degs)
                s.Degree = degs(didx);
            
                d.setModel(m);
                rmodels{didx} = d.ReducedModel;%#ok<*AGROW>

                for idx=1:m.System.InputCount
                    m.Name = sprintf('\\theta=%.4f',s.Degree);
                    d.start([],idx);
                    
                    md = d.ModelData(end);
                    % 1 = full, 7 = LSLE TD
                    absmax = 30*abs(md.ErrT(7)-md.ErrT(1));
                    relmax = 30*abs(md.RelErrT(7)-md.RelErrT(1));

                    f = d.Figures{1};
                    a = gca(f);
                    axis(a,[0 m.T md.MinErr*.9 min(1e4,absmax)]); 
                    set(legend(a),'Location','NorthEast');
                    general.Utils.saveFigure(f,sprintf('WH10_in%d_deg%f_errors',idx,s.Degree),'fig');
                    title(a,'');
                    general.Utils.saveFigure(f,sprintf('WH10_in%d_deg%f_errors',idx,s.Degree));

                    f = d.Figures{2}; a = gca(f);
                    axis(a,[0 m.T md.MinRelErr*.9 min(1,relmax)]);
                    set(legend(a),'Location','NorthEast');
                    general.Utils.saveFigure(f,sprintf('WH10_in%d_deg%f_relerr',idx,s.Degree),'fig');
                    title(a,'');
                    general.Utils.saveFigure(f,sprintf('WH10_in%d_deg%f_relerr',idx,s.Degree));

                    f = d.Figures{3}; a = gca(f);
                    general.Utils.saveFigure(f,sprintf('WH10_in%d_deg%f_ctimes',idx,s.Degree),'fig');
                    title(a,'');
                    general.Utils.saveFigure(f,sprintf('WH10_in%d_deg%f_ctimes',idx,s.Degree));
                    
                    close([d.Figures{:}]);
                end
            end
            
            save d d;
            
            range = 1:4:13;
            sort = [range (range)+1 (range)+2 (range)+3];
            d.createStatsTables(sort);
            
            %WH10Experiment.getExp1SysPlots(rmodels);
        end
        
        function getExp1SysPlots(rmodels)
            ma = ModelAnalyzer;
            ma.SingleFigures = true;
            
            % Pick 0.05 model
            r = rmodels(3);
            h = ma.analyze(r,[],3);
            title(h(2),'');
            general.Utils.saveFigure(f,'Sys_U3','eps');
        end
        
        function d = Experiment4(dim)
            
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
            s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 4;
            s.UseSVDS = dim > 10000;
            s = spacereduction.RotationDecorator(s);
            s.Degree = .01;
            s.Dims = 5;
            m.SpaceReducer = s;
            
            f = m.System.f;
            f.x0 = @(mu)ones(dim,1);
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dim,1);
            
            m.System.Inputs{1} = @(t).2*sin(t);
            %B = rnd.rand(m.dim,1);
            B = linspace(0,1,m.dim)';
            m.System.B = dscomponents.LinearInputConv(B);
            
            m.TrainingInputs = 1;
            
            m.Name = 'Sinus input, mean output mapping, larger error';
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            d = EstimatorDemo;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            d.setModel(m);
            d.start([],1);
            
            f = d.Figures{1};
            a = gca(f);
            axis(a,[0 m.T 0 10]); 
            set(legend(a),'Location','NorthEast');
            general.Utils.saveFigure(f,'WH10_e1_errors','fig');
            title(a,'');
            general.Utils.saveFigure(f,'WH10_e1_errors');
            close(f);
            
            f = d.Figures{2}; a = gca(f);
            axis(a,[0 m.T 5e-5 1]);
            set(legend(a),'Location','NorthEast');
            general.Utils.saveFigure(f,'WH10_e1_relerrors','fig');
            title(a,'');
            general.Utils.saveFigure(f,'WH10_e1_relerrors');
            close(f);
            
            f = d.Figures{3}; a = gca(f);
            general.Utils.saveFigure(f,'WH10_e1_ctimes','fig');
            title(a,'');
            general.Utils.saveFigure(f,'WH10_e1_ctimes');
            close(f);
        end
        
        function d = Experiment2(dim)
            % Experiment with POD reduction and training with all 4 inputs!
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
            f = m.System.f;
            f.x0 = @(mu)ones(dim,1);
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dim,1);
            
            %m.System.Inputs{1} = @(t).4*exp(-t/4); % alter input
            m.System.Inputs{1} = @(t).04*sin(t/3);
            m.System.Inputs{2} = @(t)exp(-abs(10-t/2)); % Gutes Beispiel!
            m.System.Inputs{3} = @(t)exp(-(10-t).^2); % Sehr interessante dynamik
            m.System.Inputs{4} = @(t)(t>10)*.1; % okay, aber nicht so spannend
            m.System.B = dscomponents.LinearInputConv(ones(m.dim,1));
            
            m.TrainingInputs = 1:4;
            
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            s = spacereduction.PODReducer;
            s.Value = 4;
            s.Mode = 'abs';
            s.UseSVDS = dim > 10000;
            m.SpaceReducer = s;
            
            d = EstimatorDemo;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            d.setModel(m);

            for idx=1:m.System.InputCount
                m.Name = 'POD4in';
                d.start([],idx);

                md = d.ModelData(end);
                % 1 = full, 7 = LSLE TD
                absmax = 30*abs(md.ErrT(7)-md.ErrT(1));
                relmax = 30*abs(md.RelErrT(7)-md.RelErrT(1));

                f = d.Figures{1};
                a = gca(f);
                axis(a,[0 m.T md.ErrT(1)*.9 min(1e4,absmax)]); 
                set(legend(a),'Location','NorthEast');
                general.Utils.saveFigure(f,sprintf('WH10_POD4_in%d_errors',idx),'fig');
                title(a,'');
                general.Utils.saveFigure(f,sprintf('WH10_POD4_in%d_errors',idx));

                f = d.Figures{2}; a = gca(f);
                axis(a,[0 m.T md.RelErrT(1)*.9 min(1,relmax)]);
                set(legend(a),'Location','NorthEast');
                general.Utils.saveFigure(f,sprintf('WH10_POD4_in%d_relerr',idx),'fig');
                title(a,'');
                general.Utils.saveFigure(f,sprintf('WH10_POD4_in%d_relerr',idx));

                f = d.Figures{3}; a = gca(f);
                general.Utils.saveFigure(f,sprintf('WH10_POD4_in%d_ctimes',idx),'fig');
                title(a,'');
                general.Utils.saveFigure(f,sprintf('WH10_POD4_in%d_ctimes',idx));

                close([d.Figures{:}]);
            end
            
            d.createStatsTables;
        end
        
        function d = Experiment3(dim)
            % Experiment with truly random component-wise functions.
            rnd = RandStream('mt19937ar','Seed',2564); %2564
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
             s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 4;
            s.UseSVDS = dim > 10000;
            s = spacereduction.RotationDecorator(s);
            s.Degree = .01;
            s.Dims = 5;
            m.SpaceReducer = s;
            
            f = m.System.f;
            f.x0 = @(mu)ones(dim,1);
            n = size(f.Centers.xi,2);
            f.Ma = -exp(-rand(dim,n));
            
            m.System.Inputs{1} = @(t).1*sin(t);
            m.System.B = dscomponents.LinearInputConv(rnd.rand(m.dim,1));
            
            this.TrainingInputs = 1;
            
            m.Name = 'Sinus input, mean output mapping, larger error';
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            d = EstimatorDemo;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            d.setModel(m);
            d.start([],1);
            
%             f = d.Figures{1};
%             a = gca(f);
%             axis(a,[0 m.T 1e-5 1e4]); 
%             set(legend(a),'Location','NorthEast');
%             general.Utils.saveFigure(f,'WH10_e3_errors','fig');
%             title(a,'');
%             general.Utils.saveFigure(f,'WH10_e3_errors');
%             close(f);
%             
%             f = d.Figures{2}; a = gca(f);
%             axis(a,[0 m.T 1e-5 1]);
%             set(legend(a),'Location','SouthWest');
%             general.Utils.saveFigure(f,'WH10_e3_relerrors','fig');
%             title(a,'');
%             general.Utils.saveFigure(f,'WH10_e3_relerrors');
%             close(f);
%             
%             f = d.Figures{3}; a = gca(f);
%             general.Utils.saveFigure(f,'WH10_e3_ctimes','fig');
%             title(a,'');
%             general.Utils.saveFigure(f,'WH10_e3_ctimes');
%             close(f);
        end
         
    end
end

