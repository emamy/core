function h = KerMorLogo
% Creates the current KerMor Logo on a new figure.
%
% @author Daniel Wirtz @date 2010-04-01
%
% @change{0,5,dw,2011-08-23} Turned this m-script into a function to be
% used i.e. in the models.ReducedModel.createImage method.
%
% This function is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
h = figure;
[X,Y] = meshgrid(-10:.5:10);
Z = exp(-(X.^2+Y.^2-2*X'*Y)/20);% + .5*exp(-(X2.^2+Y2.^2-2*X2'*Y2)/5);
surfl(X,Y,Z);
lighting gouraud;
colormap autumn;
%shading flat;
grid on;
end