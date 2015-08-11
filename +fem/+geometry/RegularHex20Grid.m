function geo = RegularHex20Grid(varargin)
% RegularHex20Grid: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-08-11
%
% @new{0,7,dw,2015-08-11} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/kermor
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing        
devperc = 0;
if length(varargin) > 3
    devperc = varargin{4};
end
g8 = fem.geometry.RegularHex8Grid(varargin{1:min(length(varargin),3)});
geo = g8.toCube20Node;

% Slightly deviate the grid
pts = geo.Nodes;
if devperc > 0
    s = RandStream('mt19937ar','Seed',1);
    pts = pts + s.rand(size(pts))*devperc;
end
geo = fem.geometry.Cube20Node(pts, geo.Elements);
end