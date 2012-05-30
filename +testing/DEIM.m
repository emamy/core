classdef DEIM
% DEIM: Tests regarding the DEIM method.
%
% See also: approx.DEIM
%
% @author Daniel Wirtz @date 2012-05-03
%
% @new{0,6,dw,2012-05-03} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        function analysis_DEIM_approx(m)
            ma = tools.ModelAnalyzer(m.buildReducedModel);
            
            mu = m.Data.ParamSamples(:,5);
%             [t, pm] = ma.compareRedFull(mu);
%             t.Format = 'tex';
%             t.saveToFile('comp_red_full.tex');
%             pm.savePlots('.',{'fig','jpg','eps'}, true);
            
            d = m.Approx;
            res = d.computeDEIMErrors(m.Data.ApproxTrainData);
            %pm = tools.PlotManager(true);
            pm = tools.PlotManager(false,1,2);
            pm.FilePrefix = 'DEIM_errors';
            d.plotDEIMErrs(res, pm);
            %pm.savePlots('.',{'fig','jpg','eps'}, true);
            pm.savePlots('.',{'fig','jpg'}, true);
            
            [etrue, EE, ED] = ma.getTrajApproxErrorDEIMEstimates(mu,[]);
%             pm = tools.PlotManager(true);
            pm = tools.PlotManager(false,1,2);
            pm.FilePrefix = 'traj_approx_err';
            ma.getTrajApproxErrorDEIMEstimates_plots(m.scaledTimes, etrue, EE, ED, pm);
%             pm.savePlots('.',{'fig','jpg','eps'}, true);
            pm.savePlots('.',{'fig','jpg'}, true);
            
            [minreq, relerrs, orders, t] = ma.getMeanRequiredErrorOrders;
            t.Format = 'tex';
            t.saveToFile('meanrequirederrororders.tex');
            
            save analysis_DEIM_approx;
        end
        
        function [t, nof, nor, o, pm] = jacobian_tests(m, pm)
            % Tis function computes the jacobians of both the full and DEIM
            % approximated system functions.
            %
            % The maximum relative error for each DEIM order and
            % jacobian at a certain snapshot vector are computed and
            % displayed.
            atd = m.Data.ApproxTrainData;
            n = min(size(atd.xi,2),500);
            sp = logical(m.System.f.JSparsityPattern);
            M = 3;
            o = round(linspace(1,m.Approx.MaxOrder,M));
            t = zeros(M,n);
            nof = t;
            nor = t;
            r = m.buildReducedModel;
            for k = 1:M
                m.Approx.Order = o(k);
                r.System.f.Order = o(k);
                for i = 1:n
                    x = atd.xi(:,i);
                    z = m.Data.W'*x;
                    x = m.Data.V*z;
                    aj = m.Approx.getStateJacobian(x,atd.ti(i),atd.mui(:,i));
                    j = m.System.f.getStateJacobian(x,atd.ti(i),atd.mui(:,i));
                    rj = r.System.f.getStateJacobian(z,atd.ti(i),atd.mui(:,i));
                    hlp = abs((aj-j)./j);
                    hlp(~sp) = 0;
                    t(k,i) = max(hlp(:));
                    nof(k,i) = norm(full(j));
                    nor(k,i) = norm(rj);
                end
            end
            if nargin < 2
                pm = tools.PlotManager(false,3,1);
            end
            h = pm.nextPlot('max_rel_err',sprintf('Maximum relative errors of DEIM jacobian vs. original one\nmasked to sparsity pattern'));
            plot(h,t')
            lbl = arrayfun(@(e)sprintf('%d',e),o,'Uniform',false);
            legend(h, lbl{:});
            
            h = pm.nextPlot('jfull_norm','Full Jacobian norms','snapshot','L2 matrix norm');
            plot(h,nof');
            legend(h, lbl{:});
            
            h = pm.nextPlot('jred_norm','Reduced Jacobian norms','snapshot','L2 matrix norm');
            plot(h,nor');
            legend(h, lbl{:});
           
            pm.done;
        end
    end    
end
