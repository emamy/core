function [Y,E] = ParamSweep(rmodel, mu, inputidx, param, paramvals)
% ParamSweep: Plots the output with error bounds for a range of one specified parameter.
%
% This function may only be used if the output of the reduced models is one-dimensional.
%
% Parameters:
% rmodel: The reduced model to use for simulations
% mu: The current parameter `\mu`. Pass this if the model has more parameters than the one to sweep.
% If set, the value at the index of the sweeped parameter will be overwritten.
% inputidx: The input index to use. As usual, set to [] if none are used.
% param: Either a string containing the name of the parameter to sweep, or it's index in the
% parameter vector.
% paramvals: A row vector of parameter values for sweeping.
%
% Return values:
% Y: The simulation matrix for every parameter, each simulation in a row.
% E: The estimated output errors for the respective setting in each row, if the reduced model's
% error estimator is enabled.
%
% @author Daniel Wirtz @date 2011-07-13
%
% @change{0,6,dw,2011-12-05} Moved this class from the \c visual/PlotParamSweep to tools.ParamSweep
%
% @new{0,5,dw,2011-07-13} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

% Validity checks 
if ischar(param)
    pidx = rmodel.getParamIndexFromName(param);
    pname = param;
elseif isposintscalar(param)
    pidx = param;
    pname = rmodel.System.Params(pidx).Name;
else
    error('param must be either a char array (parameter name) or positive parameter index');
end
npar = length(paramvals);
if npar == 0
    error('Parameter values to sweep must be given.');
end

p = 0;
fprintf('Starting parameter sweep for %s... ',pname);

mu(pidx) = paramvals(1);
[t, y] = rmodel.simulate(mu, inputidx);
if size(y,1) > 1
    error('PlotParamSweep only applicable for onedimensional system output.');
end

Y = zeros(npar,length(t));
Y(1,:) = y;
if rmodel.ErrorEstimator.Enabled
    E = zeros(npar,length(t));
    E(1,:) = rmodel.ErrorEstimator.OutputError;
else
    E = [];
end

% Iterate over parameter values
for idx = 2:npar
    mu(pidx) = paramvals(idx);
    [t, Y(idx,:)] = rmodel.simulate(mu, inputidx);
    
    if rmodel.ErrorEstimator.Enabled
        E(idx,:) = rmodel.ErrorEstimator.OutputError;
    end
    
    perc = idx/npar;
    if perc > p
        fprintf('%2.0f%% ',round(perc*100));
        p = ceil(perc*10)/10;
    end
end
fprintf('\n');

%% Prepare plot
[T,MU] = meshgrid(t,paramvals);
tit = sprintf('Outputs for parameter sweep of "%s" from %1.3f to %1.3f, base \\mu = [%s]',pname,paramvals(1),paramvals(end),num2str(mu'));

%% Upper bound
if rmodel.ErrorEstimator.Enabled
    obj = mesh(T,MU,Y+E,'EdgeColor','none','FaceColor','red');
    alpha(obj,.3);
%     mesh(T,MU,Y+E,'EdgeColor','none','FaceColor','red');
    title(tit); axis tight; xlabel('t'); ylabel(pname);
end

%% y plot
surf(T,MU,Y,'EdgeColor','none');

%% Lower bound
if rmodel.ErrorEstimator.Enabled
    figure;
    colormap jet;
    obj = surf(T,MU,Y-E,'EdgeColor','none','FaceColor','red');
    alpha(obj,.3);
%     mesh(T,MU,Y-E,'EdgeColor','none','FaceColor','red');
    hold on;
end
