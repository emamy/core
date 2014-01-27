function [pm, Y, E] = ParamSweep(rmodel, mu, inputidx, param, paramvals, pm)
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
% pm: A PlotManager instance to use for plotting
%
% Return values:
% h: The handle of the created figure
% Y: The simulation matrix for every parameter, each simulation in a row.
% E: The estimated output errors for the respective setting in each row, if the reduced model's
% error estimator is enabled.
%
% @author Daniel Wirtz @date 2011-07-13
%
% @change{0,6,dw,2011-12-05} Moved this class from the \c visual/PlotParamSweep to ParamSweep
%
% @new{0,5,dw,2011-07-13} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
if nargin < 6
    pm = PlotManager;
end

% Iterate over parameter values
Y = zeros(npar,length(rmodel.Times));
E = [];
if rmodel.ErrorEstimator.Enabled
    E = Y;
end
pi = ProcessIndicator('Starting parameter sweep with %d parameter sets for %s...',npar,false,npar,pname);
for idx = 1:npar
    mu(pidx) = paramvals(idx);
    [t, y] = rmodel.simulate(mu, inputidx);
    if size(y,1) > 1
        y = Norm.L2(y);
    end
    Y(idx,:) = y;
    
    if rmodel.ErrorEstimator.Enabled
        E(idx,:) = rmodel.ErrorEstimator.OutputError;%#ok
    end
    pi.step;
end
pi.stop;

%% Prepare plot

[T,MU] = meshgrid(t,paramvals);
ustr = '';
if ~isempty(inputidx)
    ustr = sprintf(', u_%d',inputidx);
end
tit = sprintf('Outputs for parameter sweep of "%s" (idx:%d) from %1.3f to %1.3f, base \\mu = [%s]%s',...
    pname,param,paramvals(1),paramvals(end),num2str(mu'),ustr);
h = pm.nextPlot('param_sweep',tit,'t',sprintf('Parameter %d: %s',param,pname));
%% Upper bound
if rmodel.ErrorEstimator.Enabled
    obj = surf(h,T,MU,Y+E,'EdgeColor','none','FaceColor','red');
    alpha(obj,.3);
%     mesh(T,MU,Y+E,'EdgeColor','none','FaceColor','red');
    hold on;
end

%% y plot
surf(h,T,MU,Y,'EdgeColor','none');
colormap jet;

%% Lower bound
if rmodel.ErrorEstimator.Enabled    
    obj = surf(h,T,MU,Y-E,'EdgeColor','none','FaceColor','red');
    alpha(obj,.3);
%     mesh(T,MU,Y-E,'EdgeColor','none','FaceColor','red');
    hold off;
end
if nargout < 1
    pm.done;
end