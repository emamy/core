function geo = RegularHex8Grid(xr,yr,zr,devperc)
% RegularHex8Grid: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-08-10
%
% @new{0,8,dw,2015-08-10} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/kermor
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

if nargin < 4
    devperc = 0;
    if nargin < 3
        zr = [-1 1];
        if nargin < 2
            yr = [-1 1];
            if nargin < 1
                xr = [-1 1];
            end
        end
    end
end
% Generate regular grid
[X,Y,Z] = ndgrid(xr,yr,zr);
pts = [X(:) Y(:) Z(:)]';

cubes = double.empty(0,8);
for i = 1:size(X,1)-1
    for j = 1:size(X,2)-1
        for k = 1:size(X,3)-1
            hx = X([i i+1],[j j+1],[k k+1]);
            hy = Y([i i+1],[j j+1],[k k+1]);
            hz = Z([i i+1],[j j+1],[k k+1]);
            cubes(end+1,:) = Utils.findVecInMatrix(pts,[hx(:) hy(:) hz(:)]');%#ok
        end
    end
end

% Slightly deviate the grid if wanted
if devperc > 0
    s = RandStream('mt19937ar','Seed',1);
    pts = pts + s.rand(size(pts))*devperc;
end
geo = fem.geometry.Cube8Node(pts, cubes);

end