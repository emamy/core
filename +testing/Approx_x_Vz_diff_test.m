function Approx_x_Vz_diff_test(r)
% Approx_x_Vz_diff_test(models.ReducedModel r)
%
% Computes the norm distance of `\no{\fa(x(t)) - \fa(Vz(t))} \fo t\in[0,T]`. The trajectories are
% all from the Data.TrainingData.
%
% @author Daniel Wirtz @date 2011-05-16
%
% @new{0,4,dw,2011-05-16} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
m = r.FullModel;
td = m.Data.TrainingData;
x = td(4:end,:);
t = td(3,:);
mu = td(1,:);
if any(mu ~= 0)
    mu = m.System.getParams(mu);
else
    mu = [];
end
%fax = m.Approx.evaluate(x,t,mu);
fax = m.System.f.evaluate(x,t,mu);

% Compute reduced trajectories
r.ErrorEstimator.Enabled = false;
confs = m.Data.TrainingData(1:2,1:length(m.Times):end);
mu = []; inputidx = []; xr = [];
for idx = 1:size(confs,2)
    if confs(1,idx) ~= 0
        mu = m.System.getParams(confs(1,idx));
    end
    if confs(2,idx) ~= 0
        inputidx = confs(2,idx);
    end
    [t,hlp] = r.computeTrajectory(mu,inputidx);
    xr = [xr hlp];%#ok
end
vxr = m.Data.V*xr;
favz = m.System.f.evaluate(vxr,t,mu);
%favz = m.Approx.evaluate(vxr,t,mu);

xno = sqrt(sum((x-vxr).^2));
fno = sqrt(sum((fax-favz).^2));
b = 1:size(x,2);
subplot(1,2,1);
plot(b,xno,'r',b,fno,'b');
subplot(1,2,2);
plot(b,fno./xno,'g');

end