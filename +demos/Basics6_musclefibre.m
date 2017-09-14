dim = 80;

%% Model and system setup
% The original burgers model has a nonzero initial condition and no inputs.
% we add the same things here as also been used in the a-posteriori error
% estimation paper.
m = models.musclefibre_hh.Model(dim, 500, 1, 3.828, 0.05, 5, 0.005, 1);
% This solver solves the linear part implicitly and the nonlinear party
% explicitly.
m.ODESolver = solvers.SemiImplicitEuler(m);
% This determines the plot angle for the burgers plots.
m.SaveTag = sprintf('hhTrajectory');
m.Name = sprintf('HodgkinHuxley');

s = m.System;
% Zero initial value
% Add Inputs
s.Inputs{1} = @(t)([0 0 0 0 0 0 0 0 0 2000 0 0 0 0 0 0 0 0 0 0]');
% %% Setup reduction
 m.TrainingInputs = 1;
% 
% % Sampling - log-grid
 m.Sampler = sampling.GridSampler;
% 
% % Space reduction
% % m.ComputeTrajectoryFxiData = true;
  p = spacereduction.PODReducer;
  p.IncludeBSpan = true;
  p.Mode = 'abs';
  p.Value = 40;
% % p.IncludeTrajectoryFxiData = true;
  m.SpaceReducer = p;
%             
% % Approximation of nonlinearity: Choose DEIM method here
  a = approx.DEIM(m.System);
  a.MaxOrder = 50;
  m.Approx = a;
% 
% % Run offline phase
  m.offlineGenerations;
save basic3_nonlinear;
%load basic3_nonlinear;

%% Build reduced model and do some analysis
% Details on construction of reduced models see Basics2!
r = m.buildReducedModel;
% Set the reduced DEIM approximation order to 25, and 10 orders for the
% error estimator.
r.System.f.Order = [25 10];
mu = [];
ma = ModelAnalyzer(r);
ma.compareRedFull(mu,1);
ma.analyzeError(mu,1);

%% Goodie: start the DEIM error estimator analyzer!
%DEIMEstimatorAnalyzer(r);

function [Vector] = inputOnePeak(dim, n, I_Stim, stimstart, stimend,  t)
    % t and stimstart, stimend must have the same unit of time
    nodes = dim/4;
    Vector = zeros(nodes, 1);
    if (t >= stimstart) && (t<= stimend) 
        Vector(n) = I_Stim;
    end    
end

function [Vector] = inputFrequency(dim, n, I_Stim, T, hz, abatement,  t)
    % t, 1/hz, T and abatement must have the same unit of time
    nodes = dim/4;
    Vector = zeros(nodes, 1);
    Tnew = T - abatement;
    tfreq = 0.0;
    i = 1;
    while tfreq <= Tnew
        min = tfreq;
        tfreq = tfreq + (1.0/hz);
        max = tfreq;
        tfreq = tfreq + (1.0/hz);
        freqArray(i, 1) = min;
        freqArray(i, 2) = max;
    end
    for j=1:length(freqArray)
        if (t>= freqArray(j,1)) && (t<freqArray(j,2))
            Vector(n) = I_Stim;
        end
    end
end

function [Vector] = inputFrequencyVarIStim(dim, n, I_Stim, T, hz, abatement, distribution,  t)
    % t, 1/hz, T and abatement must have the same unit of time
    nodes = dim/4;
    Vector = zeros(nodes, 1);
    Tnew = T - abatement;
    tfreq = 0.0;
    i = 1;
    while tfreq <= Tnew
        min = tfreq;
        tfreq = tfreq + (1.0/hz);
        max = tfreq;
        tfreq = tfreq + (1.0/hz);
        freqArray(i, 1) = min;
        freqArray(i, 2) = max;
    end
    for j=1:length(freqArray)
        if (t>= freqArray(j,1)) && (t<freqArray(j,2))
            Vector(n) = distribution(I_Stim);
        end
    end
end