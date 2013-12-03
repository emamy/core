function [h, Y, E] = ParamSweep2D(rmodel, mu, inputidx, params, range1, range2, ax)
% ParamSweep: Plots the output with error bounds for a range of one specified parameter.
%
% This function may only be used if the output of the reduced models is one-dimensional.
%
% Parameters:
% rmodel: The reduced model to use for simulations
% mu: The current parameter `\mu`. Pass this if the model has more parameters than the one to sweep.
% If set, the value at the index of the sweeped parameter will be overwritten.
% inputidx: The input index to use. As usual, set to [] if none are used.
% params: The indices of the parameters to sweep in the parameter vector. @type matrix<integer>
% range1: The range for parameter one @type rowvec<double>
% range2: The range for parameter two @type rowvec<double>
% ax: An axes handle to plot to. @type handle @default new figure
%
% Return values:
% h: The handle of the created figure
% Y: The simulation value at end time T over the parameter mesh.
% E: The estimated output errors for end time T over the parameter mesh.
%
% @author Daniel Wirtz @date 2011-07-13
%
% @change{0,6,dw,2011-12-05} Moved this class from the \c visual/PlotParamSweep to ParamSweep
%
% @new{0,5,dw,2011-07-13} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

% Validity checks 
[pname{1:2}] = rmodel.System.Params(params).Name;
if size(range1,2) == 0 || size(range2,2) == 0
    error('Parameter ranges to sweep must be given.');
end
if nargin < 7
    h = figure;
    ax = gca(h);
end

p = 0;
[MU1, MU2] = meshgrid(range1,range2);

pi = ProcessIndicator('Starting parameter sweep for "%s" and "%s" (%d runs)... ',numel(MU1),false,...
    pname{1},pname{2},numel(MU1));

mu(params(1)) = range1(1);
mu(params(2)) = range2(1);
[~, y] = rmodel.simulate(mu, inputidx);
if size(y,1) > 1
    error('PlotParamSweep only applicable for onedimensional system output.');
end

Y = zeros(size(MU1));
Y(1,1) = y(end);
if rmodel.ErrorEstimator.Enabled
    E = zeros(size(MU1));
    E(1,1) = rmodel.ErrorEstimator.OutputError(end);
else
    E = [];
end

% Iterate over parameter values
pi.step;
for idx = 2:numel(MU1)
    mu(params(1)) = MU1(idx);
    mu(params(2)) = MU2(idx);
    [~, y] = rmodel.simulate(mu, inputidx);
    Y(idx) = y(end);
    
    if rmodel.ErrorEstimator.Enabled
        E(idx) = rmodel.ErrorEstimator.OutputError(end);
    end
    pi.step;
end
pi.stop;

%% Prepare plot
tit = sprintf('Outputs at T=%1.2f for 2D parameter sweep, base \\mu = [%s]\n"%s" (idx:%d) from %1.3f to %1.3f\n"%s" (idx:%d) from %1.3f to %1.3f',...
    rmodel.T,num2str(mu'),pname{1},params(1),range1(1),range1(end),pname{2},params(2),range2(1),range2(end));

%% Upper bound
if rmodel.ErrorEstimator.Enabled
    obj = surf(ax,MU1,MU2,Y+E,'EdgeColor','none','FaceColor','red');
    alpha(obj,.3);
%     mesh(T,MU,Y+E,'EdgeColor','none','FaceColor','red');
    hold on;
end

%% y plot
surf(ax,MU1,MU2,Y,'EdgeColor','none');
colormap jet;
title(tit); 
xlabel(sprintf('Parameter %d: %s',params(1),pname{1}));
ylabel(sprintf('Parameter %d: %s',params(2),pname{2}));

%% Lower bound
if rmodel.ErrorEstimator.Enabled    
    obj = surf(ax,MU1,MU2,Y-E,'EdgeColor','none','FaceColor','red');
    alpha(obj,.3);
%     mesh(T,MU,Y-E,'EdgeColor','none','FaceColor','red');
    hold off;
end