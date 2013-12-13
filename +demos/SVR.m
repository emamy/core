classdef SVR
    % SVR: Support vector machine related KerMor demos
    %
    % @author Daniel Wirtz @date 2013-08-23
    %
    % @new{0,7,dw,2013-08-23} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    methods(Static)
        function pm = EpsLoss(ep)
            % Plots the epsilon-insensitive loss function for a given epsilon.
            %
            % Parameters:
            % ep: The `\epsilon` to use. @type double @default 0.1
            %
            % Return values:
            % pm: A PlotManager instance to save created figures. Optional.
            % @type PlotManager
            
            if nargin == 0
                ep = .1;
            end
            
            x = -5*ep:ep/10:5*ep;
            lx = max(0, abs(x)-ep);
            pm = PlotManager;
            h = pm.nextPlot('eps_loss',...
                sprintf('\\epsilon-insensitive loss function for \\epsilon=%g',ep),...
                'value','epsilon-norm');
            plot(h,x,zeros(size(x)),'k--',x,lx,'b',[0 0+eps],[-ep/3 max(lx)],'black--');
            
            % Create textbox
            annotation(gcf,'textbox',[0.4258 0.113 0.0461 0.07664],...
                'String',{'-\epsilon'},...
                'FontWeight','bold',...
                'FontSize',20,'LineStyle','none');
            
            
            annotation(gcf,'textbox',[0.5858 0.1229 0.0394 0.06729],...
                'String',{'\epsilon'},...
                'FontWeight','bold',...
                'FontSize',20,'LineStyle','none');
            
            axis([x(1) x(end) -ep/3 max(lx)]);
            
            if nargout < 1
                pm.done;
                pm.LeaveOpen = true;
            end
        end
        
        function ScalarEpsSVR_SMO(version, pm)
            % Demonstrates the general.regression.ScalarEpsSVR_SMO class
            %
            % Call with either 1 or 2 as argument for different SMO
            % strategies. 2 uses 2D SMO, which results in significantly
            % faster computation time with similar results.
            %
            % Parameters:
            % version: The version to use. Either 1 or 2 for 1D or 2D SMO
            % @type integer @default 2
            % pm: A PlotManager instance to use for figure creation. @type
            % PlotManager @default PlotManager
            %
            % @new{0,5,dw,2011-10-05} Added this function.
            %
            % This class is part of the framework
            % KerMor - Model Order Reduction using Kernels:
            % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
            % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
            % - \c License @ref licensing
            
            if nargin < 2
                pm = PlotManager;
                pm.LeaveOpen = true;
                if nargin < 1
                    version = 2;
                end
            end
            
            x = -5:.1:5;
%             fx = sin(x)+.2*x;
            fx = sin(x-.5).*x; 
            fx = fx ./ max(abs(fx));
            
            svr = general.regression.ScalarEpsSVR_SMO;
            svr.Version = version;
            svr.Eps = .1;
            svr.Lambda = 1/20;%1/20; % i.e. C=10 as in ScalarEpsSVR
            svr.Vis = 1;
            
            kernel = kernels.GaussKernel(.8);
            svr.K = kernel.evaluate(x,x);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernel;
            [ai, svidx] = svr.computeKernelCoefficients(fx,[]);    
            kexp.Centers.xi = x(:,svidx);
            kexp.Ma = ai';
            afx = kexp.evaluate(x);
            
            h = pm.nextPlot('svr_smo',sprintf('It: %d, #SV=%d, eps=%f',...
                svr.LastIterations,length(svidx),svr.Eps),'x','f(x)');
            demos.SVR.plotSVR(h, x, fx, afx, svr.Eps, svidx);
        end
        
        function ScalarEpsSVR_SMO_SingleSteps(varargin)
            % Demonstrates the general.regression.ScalarEpsSVR_SMO class
            %
            % Call with either Version as 1 or 2 to use different
            % SMO strategies. 2 uses 2D SMO, which results in significantly
            % faster computation time with similar results.
            %
            % Parameters:
            % varargin: Various settings for the demo.
            % Version: The version to use. Either 1 or 2 for 1D or 2D SMO
            % @type integer @default 2
            % Steps: The number of single iterations to perform. @type
            % integer @default 8
            % PM: A PlotManager instance to use for plotting. @type
            % PlotManager @default 2x2 subfigures
            % Gamma: A `\gamma` value for the kernel.
            % @type double @default 0.8
            % Eps: The `\epsilon` to use. @type double @default 0.1
            %
            % @new{0,7,dw,2013-08-23} Added this function.
            %
            % This class is part of the framework
            % KerMor - Model Order Reduction using Kernels:
            % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
            % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
            % - \c License @ref licensing
            ip = inputParser;
            ip.addParamValue('Version',2);
            ip.addParamValue('Steps',8);
            pmd = PlotManager(false,2,2);
            pmd.LeaveOpen = true;
            ip.addParamValue('PM',pmd);
            ip.addParamValue('Gamma',.8);
            ip.addParamValue('Eps',.1);
            ip.parse(varargin{:});
            r = ip.Results;
            
            x = -5:.1:5;
            fx = sin(x-.5).*x;
            fx = fx ./ max(abs(fx));
            
            svr = general.regression.ScalarEpsSVR_SMO;
            svr.Version = r.Version;
            svr.Eps = r.Eps;
            svr.Lambda = 1/20;
            svr.Vis = 1;
            
            kernel = kernels.GaussKernel(r.Gamma);
            svr.K = kernel.evaluate(x,x);
            kexp = kernels.KernelExpansion;
            kexp.Kernel = kernel;
            
            for l = 1:r.Steps
                svr.MaxCount = l;
                
                [ai, svidx] = svr.computeKernelCoefficients(fx,[]);    
                kexp.Centers.xi = x(:,svidx);
                kexp.Ma = ai';
                afx = kexp.evaluate(x);
                
                h = r.PM.nextPlot(sprintf('iter%d',l),...
                sprintf('Iteration %d',l),'x','f(x) and approximation');
                demos.SVR.plotSVR(h, x, fx, afx, svr.Eps, svidx);
            end
        end
        
    end
    
    methods(Static,Access=private)
        function plotSVR(h, x, fx, afx, eps, svidx)
            plot(h,x,fx,'r',x,[fx-eps; fx+eps],'r--');
            hold(h,'on');
            plot(h,x,afx,'b');
            % No need for eps-margin on approximation
%             plot(h,x,[afx-eps; afx+eps],'b--');
            skipped = setdiff(1:length(x),svidx);
            plot(h,x(svidx),fx(svidx),'k.','MarkerSize',16);
%             plot(h,x(skipped),fx(skipped),'xr');
        end
    end
    
end