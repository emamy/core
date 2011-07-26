function [m, r, d] = ICIAMMain
    % Current version works with KerMor 0.5

%     clear classes;
    dim = 1000;
    
    mu = [.5; 1];
    in = 1;
    
    m = ICIAMExperiment(dim);
    m.offlineGenerations;
    r = m.buildReducedModel;
    
    d = EstimatorDemo;
    d.EstimatorIterations = [1 2];
    d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
    %d.SingleFigures = true;
    d.setModel(r);
%     
    d.start(mu,in);
%     
%     ma = ModelAnalyzer;
%     ma.SingleFigures = d.SingleFigures;
%     ma.analyze(r, mu, in);
    
%     PlotParamSweep(r, mu, in, 1, 0:.1:.2); %-.9:.5:1.5
    PlotParamSweep(r, mu, in, 1, 0:.1:1);
    
    %% Experiments regardign the initial value
    PlotParamSweep(r, [1; 1], in, 2, 0.25:.1:1.2);
    PlotParamSweep(r, [.5; 1], in, 2, 0.25:.1:1.2);
    PlotParamSweep(r, [0; 1], in, 2, 0.25:.1:1.2);
    % Extranah
%     PlotParamSweep(r, [.7; 1], 1, 2, 0.9:.001:.94);
end