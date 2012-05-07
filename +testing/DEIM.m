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
    end
    
end