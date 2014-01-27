% Basics4_NonlinearKernel: 
%
% This covers:
% - Nonlinear dynamical system with nonlinearity approximation by kernel
% methods (VKOGA)
%
% @author Daniel Wirtz @date 2013-12-18
%
% @new{0,7,dw,2013-12-18} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

m = models.circ.RCLadder;
[t,y] = m.simulate([],2);
m.plot(t,y);
title('Plot for full model');

m.offlineGenerations;
r = m.buildReducedModel;
[t,yr] = r.simulate([],2);
m.plot(t,yr);
title('Plot for reduced model');
m.plot(t,yr-y);
title('Error plot for reduced/full model');

save basics4_nonlinearkernel;

%% Some analysis
ma = ModelAnalyzer(r);
ma.analyzeError([],1);

% Für advanced:
% print table
%% DPCM
%% Devel
%% Trajectory Caching
% real time plotting
