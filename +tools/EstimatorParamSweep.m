function [Y,E] = EstimatorParamSweep(ea, mu, ~, param, paramvals)
% ParamSweep: Plots the output with error bounds for a range of one specified parameter.
%
% This function may only be used if the output of the reduced models is one-dimensional.
%
% Parameters:
% ea: A complete tools.EstimatorAnalyzer instance with reduced model already set
% mu: The current parameter `\mu`. Pass this if the model has more parameters than the one to sweep.
% If set, the value at the index of the sweeped parameter will be overwritten.
% param: Either a string containing the name of the parameter to sweep, or it's index in the
% parameter vector.
% paramvals: A row vector of parameter values for sweeping. @type rowvec<double>
%
% Return values:
% Y: The simulation matrix for every parameter, each simulation in a row.
% E: The estimated output errors for the respective setting in each row, if the reduced model's
% error estimator is enabled.
%
% @author Daniel Wirtz @date 2011-07-13
%
% @new{0,5,dw,2011-07-13} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

% Validity checks
if ~isa(ea,'tools.EstimatorAnalyzer')
    error('The first argument must be a tools.EstimatorAnalyzer instance.');
end
if ischar(param)
    pidx = ea.ReducedModel.getParamIndexFromName(param);
    pname = param;
elseif isposintscalar(param)
    pidx = param;
    pname = ea.ReducedModel.System.Params(pidx).Name;
else
    error('param must be either a char array (parameter name) or positive parameter index');
end
npar = length(paramvals);
if npar == 0
    error('Parameter values to sweep must be given.');
end

% % Iterate over parameter values
% mu(pidx) = paramvals(1);
% [errs, ctimes, relerrs] = ea.getErrorEstimates(mu, inputidx, true);
% fprintf('Estimator analyzer parameter sweep 1/%d with %s=%f\n',npar,pname,mu(pidx));
% for idx = 2:npar
%     mu(pidx) = paramvals(idx);
%     [errs(:,:,idx), ctimes(idx,:), relerrs(:,:,idx)] = ea.getErrorEstimates(mu, inputidx, true);
% 
%     fprintf('Estimator analyzer parameter sweep %d/%d with %s=%f\n',idx,npar,pname,mu(pidx));
% end
load tmp;

%% Prepare plot
[T, MU] = meshgrid(ea.Model.Times,paramvals);
tit = sprintf('Error estimator analysis for parameter sweep of "%s" from %1.3f to %1.3f, base \\mu = [%s]',...
    pname,paramvals(1),paramvals(end),num2str(mu'));

% %% Upper bound
% if rmodel.ErrorEstimator.Enabled
%     obj = mesh(T,MU,Y+E,'EdgeColor','none','FaceColor','red');
%     alpha(obj,.3);
% %     mesh(T,MU,Y+E,'EdgeColor','none','FaceColor','red');
%     title(tit); axis tight; xlabel('t'); ylabel(pname);
% end

%% y plot
p = figure;
ax = gca(p);
rotate3d(ax,'on');
view(ax,9,16);
inv = figure('Visible','off');
iax = gca(inv);
hold on;
% cm = jet;
cd = cell.empty;
for idx = 1:size(errs,1)%#ok
    Z = squeeze(errs(idx,:,:))';
    %surf(ax,T,MU,Z,'EdgeColor','none','FaceColor',cm(idx*floor(size(cm,1)/size(errs,1)),:));
    h = surf(iax,T,MU,Z,'EdgeColor','none');
    cd{idx} = get(h,'CData');
    set(h,'Parent',ax);
end
title(tit);
emin = min(errs(:));
%axis(ax,[0 ea.Model.T 1 npar emin*.9 max(emin,1e4)]);
xlabel('Time t'); ylabel(sprintf('Parameter %s value',pname)); zlabel('Error estimates');
zlim(ax,[emin*.9 max(emin,1e4)]);

cm = jet;
colormap(ax,cm);
for idx = 1:size(errs,1)
    ncd = cd{idx};
    %ncd = ncd-min(ncd(:));
    ncd = size(cm,1)*max(ncd(:))./ncd;
    set(h,'CData',ncd,'CDataMapping','direct');
end
close(inv);

% %% Lower bound
% if rmodel.ErrorEstimator.Enabled
%     figure;
%     colormap jet;
%     obj = surf(T,MU,Y-E,'EdgeColor','none','FaceColor','red');
%     alpha(obj,.3);
% %     mesh(T,MU,Y-E,'EdgeColor','none','FaceColor','red');
%     hold on;
% end
